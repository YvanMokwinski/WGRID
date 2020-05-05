#include "wmesh-status.h"

extern "C"
{

    const char * wmesh_status_to_string(wmesh_status_t that_)
    {
    switch(that_)
      {
      case WMESH_STATUS_SUCCESS : return "success";
      case WMESH_STATUS_ERROR_MEMORY 	 : return "error_memory";
      case WMESH_STATUS_INVALID_ARGUMENT 	 : return "invalid_argument";
      case WMESH_STATUS_INVALID_WORK_SIZE 	 : return "invalid_work_size";
      case WMESH_STATUS_INVALID_SIZE 	 : return "invalid_size";
      case WMESH_STATUS_INVALID_POINTER 	 : return "invalid_pointer";
      case WMESH_STATUS_INVALID_CHECK 	 : return "invalid_check";
      case WMESH_STATUS_NOT_IMPLEMENTED 	 : return "not_implemented";
      case WMESH_STATUS_INVALID_CONFIG 	 : return "invalid_config";
      case WMESH_STATUS_INVALID_ENUM 	 : return "invalid_enum";
      case WMESH_STATUS_ERROR_WORKSPACE 	 : return "error_workspace";
      }
    return nullptr;
  };

};
