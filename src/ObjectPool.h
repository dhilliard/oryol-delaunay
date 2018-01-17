#pragma once

#include "Core/Containers/Array.h"

template <typename TYPE>
class ObjectPool {
public:
	inline size_t Distance(TYPE & value) const {
		return &value storage.begin();
	}
	template <typename U> inline size_t Distance(U & value) {
		return &value - reinterpret_cast<U*>(storage.begin());
	}
	template <typename ... Args>
	TYPE & Request(Args&& ... args);
	void Release(TYPE & object);
	size_t Size() const { return storage.Size() - freeSlots.Size(); }
	TYPE & operator[](size_t index) { return storage[index]; }
	const TYPE & operator[](size_t index) const { return storage[index]; }
	void Clear() {
		storage.Clear();
		freeSlots.Clear();
	}
	ObjectPool();
private:
;
	Oryol::Array<TYPE> storage;
	Oryol::Queue<size_t> freeSlots;

};
template<typename TYPE> ObjectPool<TYPE>::ObjectPool() {}

template<typename TYPE> template<typename ... Args> TYPE & ObjectPool<TYPE>::Request(Args&& ... args) {
	uint32_t index;
	if (freeSlots.Empty()) {
		index = storage.Size();
		storage.Add(args...);
	}
	else {
		index = freeSlots.Dequeue();
	}
	return storage[index];
}

template<typename TYPE> void ObjectPool<TYPE>::Release(TYPE & object) {
	freeSlots.Enqueue(Distance(object));
}