//
// Copyright (c) 2021-2024 INRIA
//

#ifndef __pinocchio_utils_openmp_hpp__
#define __pinocchio_utils_openmp_hpp__

#include <cstdlib>
#include <exception>
#include <mutex>
#include <omp.h>

namespace pinocchio
{

  /// \brief Returns the number of thread defined by the environment variable OMP_NUM_THREADS.
  ///        If this variable is not defined, this simply returns the default value 1.
  ///
  inline int getOpenMPNumThreadsEnv()
  {
    int num_threads = 1;
    
    if(const char* env_p = std::getenv("OMP_NUM_THREADS"))
      num_threads = atoi(env_p);

    return num_threads;
  }
  
  struct OpenMPException
  {
    explicit OpenMPException(const bool throw_on_deletion = false)
    : m_exception_ptr(nullptr)
    , m_throw_on_deletion(throw_on_deletion)
    {}
    
    void rethrowException() const
    {
      if(this->m_exception_ptr)
        std::rethrow_exception(this->m_exception_ptr);
    }
    
    std::exception_ptr getException() const
    {
      return m_exception_ptr;
    }
    
    bool hasThrown() const
    {
      return this->m_exception_ptr != nullptr;
    }
    
    template <typename Function, typename... Parameters>
    void run(Function f, Parameters... params)
    {
      try
      {
        f(params...);
      }
      catch (...)
      {
        this->captureException();
      }
    }
    
    void captureException() 
    {
      std::unique_lock<std::mutex> guard(this->m_lock);
      this->m_exception_ptr = std::current_exception();
    }
    
    void reset()
    {
      this->m_exception_ptr = nullptr;
    }
    
    ~OpenMPException()
    {
      if(m_throw_on_deletion)
        this->rethrowException();
    }
    
  protected:
    
    std::exception_ptr m_exception_ptr;
    std::mutex         m_lock;
    const bool         m_throw_on_deletion;
  }; // struct OpenMPException
}

#endif // ifndef __pinocchio_utils_openmp_hpp__
