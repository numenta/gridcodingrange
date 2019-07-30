/* ---------------------------------------------------------------------
 * Copyright (C) 2019, Numenta, Inc.  Unless you have an agreement
 * with Numenta, Inc., for a separate license for this software code, the
 * following terms and conditions apply:
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero Public License version 3 as
 * published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero Public License for more details.
 *
 * You should have received a copy of the GNU Affero Public License
 * along with this program.  If not, see http://www.gnu.org/licenses.
 *
 * http://numenta.org/licenses/
 * ---------------------------------------------------------------------
 */

/**
 * @file
 * Utilities for logging
 * 
 * This is a pure-header version of nupic.core's logging.
 * It has all the same functionality, except that the LogItem's output stream
 * is no longer configurable (because that would require storing it somewhere, 
 * which would require a compilation unit, and then it would not be a 
 * pure-header library.)
 */

#ifndef NTA_LOGGING_HPP
#define NTA_LOGGING_HPP

#include <iostream>
#include <sstream>

#include <stdexcept>
#include <string>
#include <utility>

#include <sstream>
#include <vector>

//----------------------------------------------------------------------

namespace ntalog
{
/**
 * @b Responsibility
 *  The Exception class is the standard Numenta exception class.
 *  It is responsible for storing rich error information:
 *  the filename and line number where the exceptional situation
 *  occured and a message that describes the exception.
 *
 * @b Rationale:
 *  It is important to store the location (filename, line number) where
 *  an exception originated, but no standard excepton class does it.
 *  The stack trace is even better and brings C++ programming to the
 *  ease of use of languages like Java and C#.
 *
 * @b Usage:
 *  This class may be used directly by instatiating an instance
 *  and throwing it, but usually you will use the NTA_THROW macro
 *  that makes it much simpler by automatically retreiving the __FILE__
 *  and __LINE__ for you and also using a wrapping LogItem that allows
 *  you to construct the exception message conveniently using the <<
 *  stream operator.
 *
 * @b Notes:
 *  1. Exception is a subclass of the standard std::runtime_error.
 *  This is useful if your code needs to interoperate with other
 *  code that is not aware of the Exception class, but understands
 *  std::runtime_error. The what() method will return the exception message
 *  and the location information will not be avaialable to such code.
 *
 *  2. Source file and line number information is useful of course
 *  only if you have access to the source code. It is not recommended
 *  to display this information to users most of the time.
 */
class Exception : public std::runtime_error
{
public:
  /**
   * Constructor
   *
   * Take the filename, line number and message
   * and store it in private members
   *
   * @param filename [const std::string &] the name of the source file
   *        where the exception originated.
   * @param lineno [unsigned] the line number in the source file where
   *        the exception originated.
   *
   * @param message [const std::string &] the description of exception
   */
  Exception(std::string filename, unsigned lineno, std::string message,
            std::string stacktrace = "")
      : std::runtime_error(""), filename_(std::move(filename)), lineno_(lineno),
        message_(std::move(message)), stackTrace_(std::move(stacktrace)) {}

  /**
   * Destructor
   *
   * Doesn't do anything, but must be present
   * because the base class std::runtime_error
   * defines a pure virtual destructor and the
   * default destructor provided by the compiler
   * doesn't have a empty exception specification
   *
   */
  virtual ~Exception() throw() {}

  /**
   * Get the exception message via what()
   *
   * Overload the what() method of std::runtime_error
   * and returns the exception description message.
   * The emptry exception specification is required because
   * it is part of the signature of std::runtime_error::what().
   *
   * @retval [const Byte *] the exception message
   */
  virtual const char *what() const throw()
  {
    try
    {
      return getMessage();
    }
    catch (...)
    {
      return "Exception caught in non-throwing Exception::what()";
    }
  }

  /**
   * Get the source filename
   *
   * Returns the full path to the source file, from which
   * the exception was thrown.
   *
   * @retval [const Byte *] the source filename
   */
  const char *getFilename() const { return filename_.c_str(); }

  /**
   * Get the line number in the source file
   *
   * Returns the (0-based) line number in the source file,
   * from which the exception was thrown.
   *
   * @retval [unsigned] the line number in the source file
   */
  unsigned getLineNumber() const { return lineno_; }

  /**
   * Get the error message
   *
   * @retval [const char *] the error message
   */
  virtual const char *getMessage() const { return message_.c_str(); }

  /**
   * Get the stack trace
   *
   * Returns the stack trace from the point the exception
   * was thrown.
   *
   * @retval [const Byte *] the stack trace
   */
  virtual const char *getStackTrace() const { return stackTrace_.c_str(); }

protected:
  std::string filename_;
  unsigned lineno_;
  std::string message_;
  std::string stackTrace_;

}; // end class Exception

/**
 * @b Description
 * A LogItem represents a single log entry. It contains a stream that
 * accumulates a log message, and its destructor calls the logger.
 *
 * A LogItem contains an internal stream
 * which is used for building up an application message using
 * << operators.
 *
 */

class LogItem
{
public:
  typedef enum
  {
    debug,
    info,
    warn,
    error
  } LogLevel;
  /**
   * Record information to be logged
   */
  LogItem(const char *filename, int line, LogLevel level)
      : filename_(filename), lineno_(line), level_(level), msg_("") {}

  /**
   * Destructor performs the logging
  */
  virtual ~LogItem()
  {
    std::ostream *out = &std::cerr;
    std::string slevel;
    switch (level_)
    {
    case debug:
      slevel = "DEBUG:";
      break;
    case warn:
      slevel = "WARN: ";
      break;
    case info:
      slevel = "INFO: ";
      out = &std::cout;
      break;
    case error:
      slevel = "ERR:";
      break;
    default:
      slevel = "Unknown: ";
      break;
    }

    (*out) << slevel << "  " << msg_.str();

    if (level_ == error)
      (*out) << " [" << filename_ << " line " << lineno_ << "]";

    (*out) << std::endl;
  }

  /*
   * Return the underlying stream object. Caller will use it to construct the
   * log message.
   */
  std::ostringstream &stream() { return msg_; }

protected:
  const char *filename_; // name of file
  int lineno_;           // line number in file
  LogLevel level_;
  std::ostringstream msg_;
};

class LoggingException : public Exception
{
public:
  LoggingException(const std::string &filename, unsigned lineno)
      : Exception(filename, lineno, std::string()), ss_(std::string()),
        lmessageValid_(false), alreadyLogged_(false) {}

  virtual ~LoggingException() throw()
  {
    if (!alreadyLogged_)
    {
      // Let LogItem do the work for us. This code is a bit complex
      // because LogItem was designed to be used from a logging macro
      LogItem *li = new LogItem(filename_.c_str(), lineno_, LogItem::error);
      li->stream() << getMessage();
      delete li;

      alreadyLogged_ = true;
    }
  }

  const char *getMessage() const override
  {
    // Make sure we use a persistent string. Otherwise the pointer may
    // become invalid.
    // If the underlying stringstream object hasn't changed, don't regenerate
    // lmessage_. This is important because if we catch this exception a second
    // call to exception.what() will trash the buffer returned by a first call
    // to exception.what()
    if (!lmessageValid_)
    {
      lmessage_ = ss_.str();
      lmessageValid_ = true;
    }
    return lmessage_.c_str();
  }

  // for Index.hpp: // because stringstream cant take << vector
  LoggingException &
  operator<<(std::vector<unsigned, std::allocator<unsigned>> v)
  {
    lmessageValid_ = false;
    ss_ << "[";
    for (auto &elem : v)
      ss_ << elem << " ";
    ss_ << "]";
    return *this;
  }

  template <typename T>
  LoggingException &operator<<(const T &obj)
  {
    // underlying stringstream changes, so let getMessage() know
    // to regenerate lmessage_
    lmessageValid_ = false;
    ss_ << obj;
    return *this;
  }

  LoggingException(const LoggingException &l)
      : Exception(l), ss_(l.ss_.str()), lmessage_(""), lmessageValid_(false),
        alreadyLogged_(true) // copied exception does not log

  {
    // make sure message string is up to date for debuggers.
    getMessage();
  }

private:
  std::stringstream ss_;
  mutable std::string lmessage_; // mutable because getMesssage() modifies it
  mutable bool lmessageValid_;
  bool alreadyLogged_;
}; // class LoggingException

} // namespace ntalog

#define NTA_DEBUG \
  ntalog::LogItem(__FILE__, __LINE__, ntalog::LogItem::debug).stream()

// Can be used in Loggable classes
#define NTA_LDEBUG(level)  \
  if (logLevel_ < (level)) \
  {                        \
  }                        \
  else                     \
    ntalog::LogItem(__FILE__, __LINE__, ntalog::LogItem::debug).stream()

// For informational messages that report status but do not indicate that
// anything is wrong
#define NTA_INFO \
  ntalog::LogItem(__FILE__, __LINE__, ntalog::LogItem::info).stream()

// For messages that indicate a recoverable error or something else that it may
// be important for the end user to know about.
#define NTA_WARN \
  ntalog::LogItem(__FILE__, __LINE__, ntalog::LogItem::warn).stream()

// To throw an exception and make sure the exception message is logged
// appropriately
#define NTA_THROW throw ntalog::LoggingException(__FILE__, __LINE__)

// The difference between CHECK and ASSERT is that ASSERT is for
// performance critical code and can be disabled in a release
// build. Both throw an exception on error.

#define NTA_CHECK(condition) \
  if (condition)             \
  {                          \
  }                          \
  else                       \
    NTA_THROW << "CHECK FAILED: \"" << #condition << "\" "

#ifdef NTA_ASSERTIONS_ON

#define NTA_ASSERT(condition) \
  if (condition)              \
  {                           \
  }                           \
  else                        \
    NTA_THROW << "ASSERTION FAILED: \"" << #condition << "\" "

#else

// NTA_ASSERT macro does nothing.
// The second line should never be executed, or even compiled, but we
// need something that is syntactically compatible with NTA_ASSERT
#define NTA_ASSERT(condition) \
  if (1)                      \
  {                           \
  }                           \
  else                        \
    ntalog::LogItem(__FILE__, __LINE__, ntalog::LogItem::debug).stream()

#endif // NTA_ASSERTIONS_ON

#endif // NTA_LOGGING_HPP
