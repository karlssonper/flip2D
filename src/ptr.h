#ifndef PTR_H_
#define PTR_H_

template <class T>
class SmartPtrInterface {
public:
    SmartPtrInterface() : ref_(0) {}
    unsigned long references() const { return ref_; }
    // DRC - support for templates
    inline const SmartPtrInterface * newRef() const { ++ref_; return this; }
    inline void deleteRef() const { if( --ref_ == 0 ) onZeroReferences(); }
protected:

    virtual ~SmartPtrInterface() {}
    virtual void onZeroReferences() const { delete this; }
private:
    mutable long unsigned ref_;
};


template <class T>
class SmartPtr
{
public:
    SmartPtr(T* p = 0) : ptr_(p) { if (ptr_) ptr_->newRef(); }
    SmartPtr(const SmartPtr<T>& mp) : ptr_(mp.ptr_) { if (ptr_) ptr_->newRef(); }
    ~SmartPtr() { if (ptr_) ptr_->deleteRef(); }

    SmartPtr<T>& operator=( const SmartPtr<T>& mp );
    SmartPtr<T>& operator=( SmartPtr<T>& mp );
    SmartPtr<T>& operator=( T* p );

    bool operator==( const SmartPtr<T>& mp ) const { return ptr_ == mp.ptr_; }
    bool operator!=( const SmartPtr<T>& mp ) const { return ptr_ != mp.ptr_; }
    bool operator==( T* p ) const { return ptr_ == p; }
    bool operator!=( T* p ) const { return ptr_ != p; }

    const T * operator->() const { return ptr_; }
    T * operator->() { return ptr_; }
    T * ptr() const { return ptr_; }

    template <class OtherType>
    operator SmartPtr<OtherType>() const { return SmartPtr<OtherType>( ptr_ ); }

    struct PointerConversion { int valid; };
    operator int PointerConversion::*() const {
        return ptr_ ? &PointerConversion::valid : 0;
    }

protected:
    T *ptr_;
};

template<class T>
SmartPtr<T>& SmartPtr<T>::operator=( const SmartPtr<T>& mp ) {
    const T * save = ptr_;
    ptr_ = mp.ptr_; 
    if( ptr_ ) ptr_->newRef();
    if( save ) save->deleteRef();
    return *this;
}

template<class T>
SmartPtr<T>& SmartPtr<T>::operator=( SmartPtr<T>& mp ) {
    T * save = ptr_;
    ptr_ = mp.ptr_; 
    if( ptr_ ) ptr_->newRef();
    if( save ) save->deleteRef();
    return *this;
}

template<class T>
SmartPtr<T>& SmartPtr<T>::operator=( T* p ) {
    T * save = ptr_;
    ptr_ = p;
    if( ptr_ ) ptr_->newRef();
    if( save ) save->deleteRef();
    return *this;
}

template <class T, class U>
SmartPtr<T> ptr_cast(SmartPtr<U> mp) {
    return dynamic_cast<T*>(mp.ptr());
}

#endif
