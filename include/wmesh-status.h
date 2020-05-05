#ifndef WMESH_STATUS_H
#define WMESH_STATUS_H

#include "wmesh-types.h"

#include <stdio.h>

//! @brief Type of the WMESH status.
typedef wmesh_int_t wmesh_status_t;

#define WMESH_STATUS_SUCCESS 		0
#define WMESH_STATUS_ERROR_MEMORY 	1
#define WMESH_STATUS_INVALID_ARGUMENT 	2
#define WMESH_STATUS_INVALID_WORK_SIZE 	3
#define WMESH_STATUS_INVALID_SIZE 	4
#define WMESH_STATUS_INVALID_POINTER 	5
#define WMESH_STATUS_INVALID_CHECK 	6
#define WMESH_STATUS_NOT_IMPLEMENTED 	7
#define WMESH_STATUS_INVALID_CONFIG 	8
#define WMESH_STATUS_INVALID_ENUM 	9
#define WMESH_STATUS_ERROR_WORKSPACE 	10

#ifdef __cplusplus
extern "C"
{
#endif
  
  const char * wmesh_status_to_string(wmesh_status_t that_);
  
#ifdef __cplusplus
}
#endif


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


#define WMESH_POINTER_CHECK(_p)  do { if (!_p) { fprintf(stderr,"WMESH invalid pointer (line=%d,file='%s')\n",__LINE__,__FILE__); return WMESH_STATUS_INVALID_POINTER; } } while(0)

#define WMESH_POSITIVE_CHECK(_p)  do { if (_p<0) { fprintf(stderr,"WMESH invalid argument (line=%d,file='%s')\n",__LINE__,__FILE__); return WMESH_STATUS_INVALID_ARGUMENT; } } while(0)

#define WMESH_STATUS_EXPECT(_expected_status,_status) do { if (_expected_status != _status) { fprintf(stderr,"WMESH failed expected status '%s', returned '%s' instead.\n",wmesh_status_to_string(_expected_status),wmesh_status_to_string(_status)); return _status; } } while(0)

#define WMESH_CHECK(_b)  (void)0 // do { if (!(_b)) { fprintf(stderr,"WMESH invalid condition (line=%d,file='%s','%s')\n",__LINE__,__FILE__,#_b); return WMESH_STATUS_INVALID_CONFIG; } } while(0)



#endif
