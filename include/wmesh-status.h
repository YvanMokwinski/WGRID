#ifndef WMESH_STATUS_H
#define WMESH_STATUS_H

#include "wmesh-types.h"

#include <stdio.h>

//! @brief Defines the status type.
typedef wmesh_int_t wmesh_status_t;

//! @brief Indicates success status.
#define WMESH_STATUS_SUCCESS 		0
//! @brief Indicates error memory status.
#define WMESH_STATUS_ERROR_MEMORY 	1
//! @brief Indicates invalid argument status.
#define WMESH_STATUS_INVALID_ARGUMENT 	2
//! @brief Indicates invalid work size status.
#define WMESH_STATUS_INVALID_WORK_SIZE 	3
//! @brief Indicates invalid size status.
#define WMESH_STATUS_INVALID_SIZE 	4
//! @brief Indicates invalid pointer status.
#define WMESH_STATUS_INVALID_POINTER 	5
//! @brief Indicates invalid check status.
#define WMESH_STATUS_INVALID_CHECK 	6
//! @brief Indicates not implemented status.
#define WMESH_STATUS_NOT_IMPLEMENTED 	7
//! @brief Indicates invalid config status.
#define WMESH_STATUS_INVALID_CONFIG 	8
//! @brief Indicates invalid enum status.
#define WMESH_STATUS_INVALID_ENUM 	9
//! @brief Indicates error workspace status.
#define WMESH_STATUS_ERROR_WORKSPACE 	10

#ifdef __cplusplus
extern "C"
{
#endif

  //! @brief Converts a status to a string.
  //! @param self_ The status.
  //! @return self_ The string of the status.
  const char * wmesh_status_to_string(wmesh_status_t self_);
  
#ifdef __cplusplus
}
#endif

//! @brief Checks a given status.
//! @param _status The status.
#define WMESH_STATUS_CHECK_EXIT(_status)				\
  if (WMESH_STATUS_SUCCESS != _status)					\
    {									\
      fprintf(stderr,"WMESH invalid status: '%s' (line=%d,file='%s')\n",wmesh_status_to_string(_status),__LINE__,__FILE__); \
      exit(_status);							\
    } (void)0


//! @brief Checks a given status.
//! @param _status The status.
#define WMESH_STATUS_CHECK(_status)					\
  do									\
    {									\
      if (WMESH_STATUS_SUCCESS != _status)				\
	{								\
	  fprintf(stderr,"WMESH invalid status: '%s' (line=%d,file='%s')\n",wmesh_status_to_string(_status),__LINE__,__FILE__);	\
	  return _status;						\
	  								\
	}								\
    } while(0)


//! @brief Checks if the pointer is not null.
//! @param _p The pointer.
#define WMESH_CHECK_POINTER(_p)  do { if (!_p) { fprintf(stderr,"WMESH invalid pointer (line=%d,file='%s')\n",__LINE__,__FILE__); return WMESH_STATUS_INVALID_POINTER; } } while(0)

//! @brief Checks if the condition is satisfied.
//! @param _c The condition.
#define WMESH_CHECK(_c)  do { if (!(_c)) { fprintf(stderr,"WMESH invalid condition (line=%d,file='%s','%s')\n",__LINE__,__FILE__,#_c); return WMESH_STATUS_INVALID_CONFIG; } } while(0)


//! @brief Checks if the condition is satisfied.
//! @param _c The condition.
#define WMESH_CHECK_POSITIVE(_c)  WMESH_CHECK( (_c) > 0 )

//! @brief Checks an expected status.
//! @param _expected_status The expected status.
//! @param _status The current status.
#define WMESH_STATUS_EXPECT(_expected_status,_status) do { if (_expected_status != _status) { fprintf(stderr,"WMESH failed expected status '%s', returned '%s' instead.\n",wmesh_status_to_string(_expected_status),wmesh_status_to_string(_status)); return _status; } } while(0)


//! @brief Checks if the condition is satisfied.
//! @param _c The condition.
#define WMESH_CHECK_STORAGE(_c)  WMESH_CHECK( WMESH_STORAGE_INTERLEAVE == (_c) || WMESH_STORAGE_BLOCK == (_c) )

#endif
