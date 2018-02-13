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
	template <typename ... Args>
	uint32_t Add(Args&& ... args);
    uint32_t Add(const TYPE & object);
	void Erase(uint32_t index);
    inline int Size() const { return activeIndices.Size(); }
    const Oryol::Set<uint32_t> & ActiveIndices() const {
        return activeIndices;
    }
    TYPE & operator[](uint32_t index);
    const TYPE & operator[](uint32_t index) const;
    template<typename U> U & GetAs(uint32_t index){
        static_assert(sizeof(TYPE) >= sizeof(U),"TYPE should be larger than U");
        static_assert(sizeof(TYPE) % sizeof(U) == 0, "sizeof(TYPE) should be a multiple of sizeof(U)");
        o_assert(IsSlotActive(index / (sizeof(TYPE)/sizeof(U))));
        return *(reinterpret_cast<U*>(storage.begin()) + index);
    }
    template<typename U> const U & GetAs(uint32_t index) const {
        static_assert(sizeof(TYPE) >= sizeof(U),"TYPE should be larger than U");
        static_assert(sizeof(TYPE) % sizeof(U) == 0, "sizeof(TYPE) should be a multiple of sizeof(U)");
        o_assert(IsSlotActive(index / (sizeof(TYPE)/sizeof(U))));
        return *(reinterpret_cast<const U*>(storage.begin()) + index);
    }
    void Clear();
    void Reserve(uint32_t amount);
    uint32_t ActiveIndexAtIndex(uint32_t index);
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

template<typename TYPE> template<typename ... Args> uint32_t ObjectPool<TYPE>::Add(Args&& ... args) {
	uint32_t index;
	if (freeSlots.Empty()) {
		index = storage.Size();
        storage.Add(std::forward<Args>(args)...);
        generation.Add(0);
	}
	else {
		index = freeSlots.Dequeue();
		storage[index] = TYPE(std::forward<Args>(args)...);
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
    activeIndices.Erase(index);
    freeSlots.Enqueue(index);
}

template<typename TYPE> uint32_t ObjectPool<TYPE>::ActiveIndexAtIndex(uint32_t index) {
    return activeIndices.ValueAtIndex(index);
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
