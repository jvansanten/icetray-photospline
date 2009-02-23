/**
 *@file
 *@brief Stand alone kogging system for PSInterface
 *
 * The logging system is used to print out various information
 * while running PSI. It is designed to be easy to replace
 * with the standard icecube logging system. To active the system
 * define ENABLE_PSI_LOGGING, if not defined the standard icecube 
 * system will be used. As opposed to an advanced logging system 
 * the PSI logging doesn't print filenames or line  numbers and 
 * is less configurable.
 *
 * The difference between logging and just normal printing is that
 * the level of severity (loggginglevel) of the printout is defined 
 * and that a custom message is prefixed to the message 
 * (PSI logginglevel). The different levels of severity are            
 * - trace   Simple printouts to follow the code.  (ie I'm here)
 * - debug   Debug statement (ie this is what happens)
 * - info    Information messages (ie you should know this)
 * - warn    Warning messages (ie something maybe wrong)
 * - error   Error messages (ie something is wrong)
 * - fatal   Severe errors that requires PSI to exit
 *
 * The printout is done by calling log_logginglevel( format ) where 
 * format is the regular input to printf. The buffer that the log is
 * printed to can be changed by defining PSI_LOGGING_BUFFER
 * (the default is stderr).
 *
 * A logging threshold can be set that removes printout below a chosen level 
 * by defining one of the flags: PSI_LOGGING_LEVEL_DEBUG, 
 * PSI_LOGGING_LEVEL_INFO, PSI_LOGGING_LEVEL_WARN, 
 * PSI_LOGGING_LEVEL_ERROR, PSI_LOGGING_LEVEL_FATAL.
 *
 * If no level specified DEFAULT will be chosen
 *
 * (All flags can be set by appending -Dflag to CFLAGS before compiling)
 * 
 *@author Thomas Burgess
 * 
 * (c) the IceCube Collaboration
 *
 * $Revision: 1.2 $
 * $Author$
 * $Date$
 * $Id$
 */

#ifndef __PSI_Logging_h__
#define __PSI_Logging_h__

#ifdef ENABLE_PSI_LOGGING

#include <cstdio>

//If no buffer specified use stderr
#ifndef PSI_LOGGING_BUFFER
#define PSI_LOGGING_BUFFER stderr
#endif

// If level specified, define all levels above it
// If no level specified level is set to default (which is DEBUG)

#ifdef PSI_LOGGING_LEVEL_TRACE
#define PSI_LOGGING_LEVEL_DEBUG
#endif

#ifdef PSI_LOGGING_LEVEL_DEBUG
#define PSI_LOGGING_LEVEL_INFO
#endif

#ifdef PSI_LOGGING_LEVEL_INFO
#define PSI_LOGGING_LEVEL_WARN
#endif

#ifdef PSI_LOGGING_LEVEL_WARN
#define PSI_LOGGING_LEVEL_ERROR
#endif

#ifdef PSI_LOGGING_LEVEL_ERROR
#define PSI_LOGGING_LEVEL_FATAL
#endif

#ifdef PSI_LOGGING_LEVEL_FATAL
#else
//No level specified set default settings (DEBUG)
#define PSI_LOGGING_LEVEL_DEBUG
#define PSI_LOGGING_LEVEL_INFO
#define PSI_LOGGING_LEVEL_WARN
#define PSI_LOGGING_LEVEL_ERROR
#define PSI_LOGGING_LEVEL_FATAL
#endif

#define SET_LOGGER( name )

#ifdef PSI_LOGGING_LEVEL_TRACE
#define log_trace(format, ...)\
    fprintf( PSI_LOGGING_BUFFER, "PSI TRACE: ");\
    fprintf( PSI_LOGGING_BUFFER, format, ##__VA_ARGS__);\
    fprintf( PSI_LOGGING_BUFFER, "\n");
#else
#define log_trace(format, ...);
#endif

#ifdef PSI_LOGGING_LEVEL_DEBUG
#define log_debug(format, ...)\
    fprintf( PSI_LOGGING_BUFFER, "PSI DEBUG: ");\
    fprintf( PSI_LOGGING_BUFFER, format, ##__VA_ARGS__);\
    fprintf( PSI_LOGGING_BUFFER, "\n");
#else
#define log_debug(format, ...);
#endif

#ifdef PSI_LOGGING_LEVEL_INFO
#define log_info(format, ...)\
    fprintf( PSI_LOGGING_BUFFER, "PSI INFO: ");\
    fprintf( PSI_LOGGING_BUFFER, format, ##__VA_ARGS__);\
    fprintf( PSI_LOGGING_BUFFER, "\n");
#else
#define log_info(format, ...);
#endif

#ifdef PSI_LOGGING_LEVEL_WARN
#define log_warn(format, ...)\
    fprintf( PSI_LOGGING_BUFFER, "PSI WARN: ");\
    fprintf( PSI_LOGGING_BUFFER, format, ##__VA_ARGS__);\
    fprintf( PSI_LOGGING_BUFFER, "\n");
#else
#define log_warn(format, ...);
#endif

#ifdef PSI_LOGGING_LEVEL_ERROR
#define log_error(format, ...)\
    fprintf( PSI_LOGGING_BUFFER, "PSI ERROR: ");\
    fprintf( PSI_LOGGING_BUFFER, format, ##__VA_ARGS__);\
    fprintf( PSI_LOGGING_BUFFER, "\n");
#else
#define log_error(format, ...);
#endif

#ifdef PSI_LOGGING_LEVEL_FATAL
#define log_fatal(format, ...)\
    fprintf( PSI_LOGGING_BUFFER, "PSI FATAL: ");\
    fprintf( PSI_LOGGING_BUFFER, format, ##__VA_ARGS__);\
    fprintf( PSI_LOGGING_BUFFER, "\n");\
    exit(-1); 
#else
#define log_error(format, ...);
#endif

#else
// Or use logging instead
#include "icetray/I3Logging.h"
#endif

#endif
