// $Id: libmesh_logging.h,v 1.3 2004-02-27 22:43:47 spetersen Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#ifndef __libmesh_logging_h__
#define __libmesh_logging_h__


// The library configuration options
#include "libmesh_common.h"

// Macros for performance logging.  This allows us
// to add performance monitors to the code without
// impacting performance when performance logging
// is disabled.
#ifdef ENABLE_PERFORMANCE_LOGGING

// Note the log is in libMesh, so we need to include it.
#  include "libmesh.h"
#  define START_LOG(a,b)   { libMesh::perflog.start_event(a,b); }
#  define STOP_LOG(a,b)    { libMesh::perflog.stop_event(a,b); }
#  define PAUSE_LOG(a,b)   { libMesh::perflog.pause_event(a,b); }
#  define RESTART_LOG(a,b) { libMesh::perflog.restart_event(a,b); }

#else

#  define START_LOG(a,b)   {}
#  define STOP_LOG(a,b)    {}
#  define PAUSE_LOG(a,b)   {}
#  define RESTART_LOG(a,b) {}

#endif





#endif // #define __libmesh_logging_h__
