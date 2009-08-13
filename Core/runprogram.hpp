/*
  This file is part of the FIRE -- Flexible Image Retrieval System

  FIRE is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  FIRE is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with FIRE; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef __runprogram__hpp__
#define __runprogram__hpp__
#include <string>

/** Basic operations to run external programs */

/** Run an external program using execvp. That is, the current process
    is replaced by the program started. */
void run(::std::string commandline); 

/** Run an external program using execvp in a new thread. That is, a
    new thread is forked and the program takes this thread. */
void runInBackground(::std::string  commandline);

/** wait for all children to terminate */
void waitForAllChildren();
#endif
