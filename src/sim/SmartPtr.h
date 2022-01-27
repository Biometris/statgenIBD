#ifndef SMART_POINTER_HEADER
#define SMART_POINTER_HEADER

#include "ibdexcept.h" 

namespace ibd
{

template<class T> 
class SmartPtr 
{
public:
	SmartPtr(): ptr_cntr(new size_t(1)), ptr_data(0) { }
	SmartPtr(T* pt): ptr_cntr(new size_t(1)), ptr_data(pt) { }
	SmartPtr(const SmartPtr& h): ptr_cntr(h.ptr_cntr), ptr_data(h.ptr_data) { ++*ptr_cntr; }
	~SmartPtr();                        

	SmartPtr& operator=(const SmartPtr&); 
	operator bool() const { return ptr_data ? true : false; }
	T& operator*() const;         
	T* operator->() const;        
	T* get() const { return ptr_data; }

private:
	T*	    ptr_data;  // pointer to the data
	size_t* ptr_cntr;  // pointer to the reference counter
};

template<class T>
T& SmartPtr<T>::operator*() const 
{ 
	if (ptr_data) 
		return *ptr_data; 
	throw ibd_error("unbound SmartPtr"); 
}

template<class T>
T* SmartPtr<T>::operator->() const 
{	
	if (ptr_data) 
		return ptr_data; 
	throw ibd_error("unbound SmartPtr"); 
}

template<class T>
SmartPtr<T>& SmartPtr<T>::operator=(const SmartPtr& rhs)
{
	++*rhs.ptr_cntr;
	if (--*ptr_cntr == 0) 
	{
		delete ptr_cntr;
        delete ptr_data;
    }
    ptr_cntr = rhs.ptr_cntr;
    ptr_data = rhs.ptr_data;
    return *this;
}

template<class T> 
SmartPtr<T>::~SmartPtr()
{
	if (--*ptr_cntr == 0) 
	{
		delete ptr_cntr;
		delete ptr_data;
    }
}

template<class T> 
T* clone(const T* tp)
{
	return tp->clone();
}

}

#endif
