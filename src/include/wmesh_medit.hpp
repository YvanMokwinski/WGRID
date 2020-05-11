#pragma once

#include "wmesh.hpp"

wmesh_int_t 	wmesh_read_medit	(wmesh_t**self_,
					 bool is_binary_,
					 const char * filename,
					 ...);

wmesh_status_t 	wmesh_write_medit	(const wmesh_t* 	self_,
					 bool  			is_binary_,
					 const char * 		filename_,
					 ...);
