#pragma once

#include "wmesh.hpp"
extern "C"
{
  wmesh_int_t wmesh_read_medit(wmesh_t**self_,const char * filename,...);
}
