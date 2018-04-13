/*
Copyright 2013-2018 Denis Hilliard <denis.z.hilliard@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of this
software and associated documentation files(the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
to whom the Software is furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#pragma once

#include "Core/Containers/Array.h"
#include "Core/Containers/Queue.h"
#include "Core/Containers/Set.h"

template <typename TYPE>
class ObjectPool {
public:
    uint32_t Distance(const TYPE & o) const {
        return &o - storage.begin();
    }
    template<typename U>
    uint32_t Distance(const U & o) const {
        return &o - reinterpret_cast<const U*>(storage.begin());
    }
    uint32_t Add(const TYPE & object);
	void Erase(uint32_t index);
    inline int Size() const { return activeIndices.Size(); }
    const Oryol::Set<uint32_t> & ActiveIndices() const {
        return activeIndices;
    }
    TYPE & operator[](uint32_t index);
    const TYPE & operator[](uint32_t index) const;
    template<typename U> U & GetAs(uint32_t index) const {
        static_assert(sizeof(TYPE) >= sizeof(U),"TYPE should be larger than U");
        static_assert(sizeof(TYPE) % sizeof(U) == 0, "sizeof(TYPE) should be a multiple of sizeof(U)");
        o_assert(IsSlotActive(index / (sizeof(TYPE)/sizeof(U))));
        return *(reinterpret_cast<U*>(storage.begin()) + index);
    }
    void Clear();
    void Reserve(uint32_t amount);
    uint32_t ActiveIndexAtIndex(uint32_t index) const;
    inline bool IsSlotActive(uint32_t index) const{
        uint32_t which = occupancy[index / 32];
        return (which & (1 << (index & 31))) > 0;
    }
    inline uint32_t SlotGeneration(const uint32_t index) const {
        return generation[index];
    }

	ObjectPool();
private:
	Oryol::Array<TYPE> storage;
    Oryol::Set<uint32_t> activeIndices;
	Oryol::Queue<uint32_t> freeSlots;
    Oryol::Array<uint32_t> occupancy;
    Oryol::Array<uint32_t> generation;
    inline void enable(uint32_t index){
        uint32_t & which = occupancy[index / 32];
        which |= (1 << (index & 31));
    }
    inline void disable(uint32_t index){
        uint32_t & which = occupancy[index / 32];
        which &= ~(1 << (index & 31));
    }

};
template<typename TYPE> ObjectPool<TYPE>::ObjectPool() {}

template<typename TYPE> uint32_t ObjectPool<TYPE>::Add(const TYPE & object){
    uint32_t index = -1;
    if (freeSlots.Empty()) {
        index = storage.Size();
        storage.Add(object);
        generation.Add(0);
    }
    else {
        index = freeSlots.Dequeue();
        storage[index] = object;
        generation[index]++;
    }
    if(occupancy.Size() <= int(index / 32))
        occupancy.Add(0);
    enable(index);
    activeIndices.Add(index);
    return index;
}

template<typename TYPE> void ObjectPool<TYPE>::Erase(uint32_t index) {
    o_assert_dbg(IsSlotActive(index));
    disable(index);
    storage[index].~TYPE();
    activeIndices.Erase(index);
    freeSlots.Enqueue(index);
}

template<typename TYPE> uint32_t ObjectPool<TYPE>::ActiveIndexAtIndex(uint32_t index) const {
    //Workaround for Oryol::Set
    return const_cast<Oryol::Set<uint32_t>&>(activeIndices).ValueAtIndex(index);
}

template<typename TYPE> void ObjectPool<TYPE>::Clear() {
    activeIndices.Clear();
    freeSlots.Clear();
    storage.Clear();
    occupancy.Clear();
}
template<typename TYPE> void ObjectPool<TYPE>::Reserve(uint32_t amount){
    if(amount > (uint32_t)freeSlots.Size())
        storage.Reserve(amount - freeSlots.Size());
}

template<typename TYPE> TYPE & ObjectPool<TYPE>::operator[](uint32_t index){
    o_assert_dbg(IsSlotActive(index));
    return storage[index];
}

template<typename TYPE> const TYPE & ObjectPool<TYPE>::operator[](uint32_t index) const{
    o_assert_dbg(IsSlotActive(index));
    return storage[index];
}
