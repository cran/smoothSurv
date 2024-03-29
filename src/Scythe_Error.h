/* Scythe_Error.h
 *
 * This header provedes class definitions and implementations for
 * errors thrown by Scythe libraries.
 *
 * Scythe C++ Library
 * Copyright (C) Kevin M. Quinn, Andrew D. Martin,
 * and Daniel B. Pemstein
 *
 * This code written by:
 *
 * Kevin Quinn
 * Assistant Professor
 * Dept. of Political Science and
 * Center for Statistics and Social Sciences
 * Box 354322
 * University of Washington
 * Seattle, WA 98195-4322
 * quinn@stat.washington.edu
 *
 * Andrew D. Martin
 * Assistant Professor
 * Dept. of Political Science
 * Campus Box 1063
 * Washington University
 * St. Louis, MO 63130
 * admartin@artsci.wustl.edu
 * 
 * Daniel B. Pemstein
 * dbpemste@artsci.wustl.edu
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */

#ifndef SCYTHE_ERROR_H
#define SCYTHE_ERROR_H

//#include <exception>    /* commented on 04/06/2013, it is nowhere needed */
#include <string>
#include <sstream>
#include <iostream>
//#include <cstdlib>      /* added by AK on 26/06/2008 to provide abort() in g++ 4.3 compiler     */
                          /* again commented on 06/01/2012 after replacement of abort by R error  */

#ifndef R_NO_REMAP        /* Added on 04/06/2013 to prevent <R_ext/Error.h> included by <R.h>        */
#define R_NO_REMAP        /* from definition of error as Rf_error which causes problems on some      */
#endif                    /* compilers, see explanation on bottom of this file close to the place    */ 
                          /* where error() is commented.                                             */
#include <R.h>            /* added by AK on 06/01/2012 to provide error()                            */

/***** The following piece of code has been motivated by /usr/include/assert.h (in Debian Linux) *****/
/***** and has been added by Arnost Komarek on 14/08/2007                                        *****/
/***** ========================================================================================= *****/

/* Version 2.4 and later of GCC define a magical variable `__PRETTY_FUNCTION__'
   which contains the name of the function currently being defined.
   This is broken in G++ before version 2.6.
   C9x has a similar variable called __func__, but prefer the GCC one since
   it demangles C++ function names.  */
//# if defined __cplusplus ? __GNUC_PREREQ (2, 6) : __GNUC_PREREQ (2, 4)
# if defined __GNUC__ && __GNUC__ >= 2
#   define __AK_PRETTY_FUNCTION__ __PRETTY_FUNCTION__
# else
#  if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#   define __AK_PRETTY_FUNCTION__ __func__
#  else
#   define __AK_PRETTY_FUNCTION__ ((__const char *) 0)
#  endif
# endif

/***** ======================================================================================== *****/

//extern std::ostringstream os_AK;       /**** added by AK on 10/10/2022, global variable being defined in smoothSurvReg84.cpp  ****/

namespace SCYTHE {

        /**** This file-local variable holds the output of the last
         * scythe_exception constructed.
         ****/
        namespace {
            std::string serr;
        }

        /**** A replacement for the default terminate handler.  This outputs
         * the string held in serr before calling abort, thereby notifying
         * the user of why the program crashed.
         ****/
        inline void scythe_terminate();
  
        /**** The scythe exception abstract base class ****/
        class scythe_exception : public std::exception
        {
		public:
			scythe_exception (const std::string &head,
					  const std::string &file,
					  const std::string &function,
					  const unsigned int &line,
					  const std::string &message = "",
					  const bool &halt = false) throw ()
				:	exception(),
					head_ (head),
					file_ (file),
					function_ (function),
					line_ (line),
					message_ (message)
			{
				std::ostringstream os;
				os << head_ << " in " << file_ << ", " << function_ << ", "
					<< line_ << ": " << message_ << "!";
				serr = os.str();
				std::set_terminate(scythe_terminate);
				if (halt)
				  //std::terminate();
				  //error("std::terminate\n", "Scythe_Error.h");         /* this also causes some problems, see comments at bottom of this file*/
 				  REprintf("ERROR in SCYTHE: %s\n\n", (char*)(&serr));
			}

			scythe_exception (const scythe_exception &e) throw()
				:	exception(),
					head_ (e.head_),
					file_ (e.file_),
					function_ (e.function_),
					line_ (e.line_),
					message_ (e.message_)
			{
			}

			scythe_exception &operator= (const scythe_exception &e) throw()
			{
				head_ = e.head_;
				file_ = e.file_;
				function_ = e.function_;
				line_  = e.line_;
				message_ = e.message_;

				return *this;
			}

			virtual ~scythe_exception() throw()
			{
			}

	                //virtual const char * what() const throw()      // the whole function commented by AK on 10/10/2022
	                //{
			        //std::ostringstream os;
				//os << head_ << " in " << file_ << ", " << function_ << ", "
				//	<< line_ << ": " << message_ << "!";
				//return os.str().c_str();

                                /*** Replacement of the above code by AK on 10/10/2022 ***/
                        //      os_AK.str(""); 
			// 	os_AK << head_ << " in " << file_ << ", " << function_ << ", "
			//		<< line_ << ": " << message_ << "!";
			//	return os_AK.str().c_str();			  
			//}

		private:
			std::string head_;
			std::string file_;
			std::string function_;
			unsigned int line_;
			std::string message_;
	};


	/**** Exception class types, added as needed ****/
	class scythe_alloc_error : public scythe_exception
	{
		public:
			scythe_alloc_error (const std::string &file,
					    const std::string &function,
					    const unsigned int &line,
					    const std::string &message = "",
					    const bool &halt = false) throw()
				:	scythe_exception("SCYTHE_ALLOCATION_ERROR", file, function, line, message, halt)
			{}
	};
	
	class scythe_invalid_arg : public scythe_exception
	{
		public:
			scythe_invalid_arg (	const std::string &file,
						const std::string &function,
						const unsigned int &line,
						const std::string &message = "",
						const bool &halt = false) throw()
				:	scythe_exception("SCYTHE_INVALID ARGUMENT", file, function, line, message, halt)
			{}
	};
	
	class scythe_file_error : public scythe_exception
	{
		public:
			scythe_file_error (	const std::string &file,
						const std::string &function,
						const unsigned int &line,
						const std::string &message = "",
						const bool &halt = false) throw()
				:	scythe_exception("SCYTHE FILE ERROR", file, function, line, message, halt)
			{}
	};
	
	class scythe_conformation_error : public scythe_exception
	{
		public:
			scythe_conformation_error (	const std::string &file,
							const std::string &function,
							const unsigned int &line,
							const std::string &message = "",
							const bool &halt = false) throw()
				:	scythe_exception("SCYTHE CONFORMATION ERROR", file, function, line, message, halt)
			{}
	};
	
	class scythe_dimension_error : public scythe_exception
	{
		public:
			scythe_dimension_error (const std::string &file,
						const std::string &function,
						const unsigned int &line,
						const std::string &message = "",
						const bool &halt = false) throw()
				:	scythe_exception("SCYTHE DIMENSION ERROR", file, function, line, message, halt)
			{}
	};
	
	class scythe_null_error : public scythe_exception
	{
		public:
			scythe_null_error (	const std::string &file,
						const std::string &function,
						const unsigned int &line,
						const std::string &message = "",
						const bool &halt = false) throw()
				:	scythe_exception("SCYTHE NULL ERROR", file, function, line, message, halt)
			{}
	};
	
	class scythe_type_error : public scythe_exception
	{
		public:
			scythe_type_error (	const std::string &file,
						const std::string &function,
						const unsigned int &line,
						const std::string &message = "",
						const bool &halt = false) throw()
				:	scythe_exception("SCYTHE TYPE ERROR", file, function, line, message, halt)
			{}
	};
	
	class scythe_out_of_range_error : public scythe_exception
        {
             public:
                 scythe_out_of_range_error (const std::string &file,
                                            const std::string &function,
                                            const unsigned int &line,
                                            const std::string &message = "",
                                            const bool &halt = false) throw()
                    :   scythe_exception("SCYTHE OUT OF RANGE ERROR", file, function, line, message, halt)
                 {}
        };
	
	class scythe_convergence_error : public scythe_exception
	{
		public:
			scythe_convergence_error (	const std::string &file,
               						const std::string &function,
		                			const unsigned int &line,
				            		const std::string &message = "",
						        const bool &halt = false) throw()
				:	scythe_exception("SCYTHE CONVERGENCE ERROR", file, function, line, message, halt)
			{}
	};
	
	class scythe_range_error: public scythe_exception
	{
		public:
			scythe_range_error(	const std::string &file,
						const std::string &function,
             					const unsigned int &line,
		             			const std::string &message = "",
				        	const bool &halt = false) throw()
				:	scythe_exception("SCYTHE RANGE ERROR", file, function, line, message, halt)
			{}
	};

	class scythe_precision_error: public scythe_exception
	{
		public:
			scythe_precision_error(	const std::string &file,
 					        const std::string &function,
						const unsigned int &line,
						const std::string &message = "",
						const bool &halt = false) throw()
				:	scythe_exception("SCYTHE PRECISION ERROR", file, function, line, message, halt)
			{}
	};

	// The definition of our terminate handler described above
	inline void scythe_terminate ()
	{
    	        // VERSION 1:
	        //std::cerr << serr << std::endl;      
	        //std::cerr << std::endl;
	        //abort();

                // VERSION 2:
                //error("%s\n\n", (char*)(&serr));     // added on 20120106 to replace previsous three rows
		/* Discovered later on: On some compilers the following error occurs: 'Rf_error' is not a member of 'std::codecvt:base'. */
                /* First reported by Brian Ripley on 19/03/2013 (compiler on Mac), later found also on Karlov cluster                    */
                /* with the Redhat Linux 64 bit compiler.                                                                                */
                /* It is caused by the fact that <R_ext/Error.h> included by <R.h> defines error as Rf_error.                            */
                /* Final decission from 20130604:  Just print the message here but otherwise do nothing.                                 */
                /* In most cases, this should not cause any problems...                                                                  */
             
                // VERSION 3 (from 20130604):
  	        REprintf("ERROR in SCYTHE: %s\n\n", (char*)(&serr));
	}

} // end namspace SCYTHE

#endif /* SCYTHE_ERROR_H */
