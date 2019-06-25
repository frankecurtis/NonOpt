// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTREPORTER_HPP__
#define __NONOPTREPORTER_HPP__

#include <cstdarg>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "NonOptEnumerations.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class Report;
class FileReport;
class StreamReport;

/**
  * Reporter class
  */
class Reporter
{

 public:
  /** @name Constructors */
  //@{
  /**
    * Declare Reporter
    */
  Reporter();
  //@}

  /** @name Destructor */
  //@{
  /**
    * Delete Reporter
    */
  virtual ~Reporter();
  //@}

  /** @name Print methods */
  //@{
  /**
    * printf
    * \param[in] type is ReportType at which to print
    * \param[in] level is ReportLevel at which to print
    * \param[in] format is formatting string
    */
  virtual void printf(ReportType type,
                      ReportLevel level,
                      const char* format,
                      ...) const;
  /**
    * Print list
    * \param[in] type is ReportType at which to print
    * \param[in] level is ReportLevel at which to print
    * \param[in] format is formatting string
    * \param[in] lst is list of strings to print
    */
  virtual void printList(ReportType type,
                         ReportLevel level,
                         const char* format,
                         va_list lst) const;
  //@}

  /** @name Add methods */
  //@{
  /**
    * Add Report
    * \param[in] report is pointer to Report to add
    * \return indicator of success (true) or failure (false)
    */
  virtual bool addReport(const std::shared_ptr<Report> report);
  /**
    * Add FileReport
    * \param[in] report_name is name of FileReport to add
    * \param[in] file_name is name of file for FileReport to add
    * \param[in] type is ReportType to set for FileReport
    * \param[in] default_level is ReportLevel to set for FileReport
    * \return indicator of success (true) or failure (false)
    */
  virtual bool addFileReport(std::string report_name,
                             std::string file_name,
                             ReportType type,
                             ReportLevel level);
  //@}

  /** @name Get methods */
  //@{
  /**
    * Get Report
    * \param[in] name is name of report to get
    * \return pointer to Report with report_name, if it exists; else, it is null
    */
  virtual std::shared_ptr<Report> report(std::string name);
  //@}

  /** @name Acceptance method */
  //@{
  /**
    * Checks if (type,level) pair is accepted by some report
    * \param[in] type is ReportType of query
    * \param[in] level is ReportLevel of query
    * \return indicator of acceptable (true) or not (false)
    */
  virtual bool isAccepted(ReportType type,
                          ReportLevel level) const;
  //@}

  /** @name Flush buffer method */
  //@{
  /**
    * Flush buffer
    */
  virtual void flushBuffer() const;
  //@}

  /** @name Delete method */
  //@{
  /**
    * Delete Reports
    */
  virtual void deleteReports();
  //@}

 private:
  /** @name Default compiler generated methods
    * (Hidden to avoid implicit creation/calling.)
    */
  //@{
  /**
    * Copy constructor
    */
  Reporter(const Reporter&);
  /**
    * Overloaded equals operator
    */
  void operator=(const Reporter&);
  //@}

  /** @name Private members */
  //@{
  std::vector<std::shared_ptr<Report> > reports_; /**< List (vector) of reports */
  //@}

};  // end Reporter

/**
  * Report class
  */
class Report
{

 public:
  /** @name Constructors */
  //@{
  /**
    * Construct Report
    * \param[in] name is name of report
    * \param[in] type is ReportType of report
    * \param[in] level is ReportLevel of report
    */
  Report(std::string name,
         ReportType type,
         ReportLevel level);
  //@}

  /** @name Destructor */
  //@{
  /**
    * Delete Report
    */
  virtual ~Report();
  //@}

  /** @name Get methods */
  //@{
  /**
    * Get name
    * \return is name of Report as string
    */
  virtual std::string name();
  //@}

  /** @name Set methods */
  //@{
  /**
    * Set type and level
    * \param[in] type is ReportType of print level to set
    * \param[in] level is ReportLevel of print level to set
    */
  virtual void setTypeAndLevel(ReportType type,
                               ReportLevel level);
  //@}

  /** @name Acceptance method */
  //@{
  /**
    * Checks if (type,level) pair is accepted by report
    * \param[in] type is ReportType of query
    * \param[in] level is ReportLevel of query
    * \return indicator of success (true) or failure (false)
    */
  virtual bool isAccepted(ReportType type,
                          ReportLevel level) const;
  //@}

  /** @name Print methods */
  //@{
  /**
    * Print list
    * \param[in] type is ReportType at which to print
    * \param[in] level is ReportLevel at which to print
    * \param[in] format is formatting string
    * \param[in] lst is list of strings to print
    */
  virtual void printList(ReportType type,
                         ReportLevel level,
                         const char* format,
                         va_list lst)
  {
    printListImplementation(type, level, format, lst);
  }
  //@}

  /** @name Flush buffer method */
  //@{
  /**
    * Flush buffer
    */
  virtual void flushBuffer() { flushBufferImplementation(); }
  //@}

  /** @name Close report method */
  //@{
  /**
    * Close report
    */
  virtual void close() { closeImplementation(); }
  //@}

 protected:
  /** @name Implementation methods */
  //@{
  /**
    * Print list implementation
    * \param[in] type is ReportType at which to print
    * \param[in] level is ReportLevel at which to print
    * \param[in] format is formatting string
    * \param[in] lst is list of strings to print
    */
  virtual void printListImplementation(ReportType type,
                                       ReportLevel level,
                                       const char* format,
                                       va_list lst) = 0;
  /**
    * Flush buffer implementation
    */
  virtual void flushBufferImplementation() = 0;
  /**
    * Close implementation
    */
  virtual void closeImplementation() = 0;
  //@}

 private:
  /** @name Default compiler generated methods
    * (Hidden to avoid implicit creation/calling.)
    */
  //@{
  /**
    * Constructor with no arguments
    */
  Report();
  /**
    * Copy constructor
    */
  Report(const Report&);
  /**
    * Overloaded equals operator
    */
  void operator=(const Report&);
  //@}

  /**
   * Private members
   */
  //@{
  std::string name_;
  ReportType type_;
  ReportLevel level_;
  //@}

};  // end Report

/**
  * FileReport class
  */
class FileReport : public Report
{

 public:
  /** @name Constructors */
  //@{
  /**
    * Construct FileReport
    * \param[in] name is name of report
    * \param[in] type is ReportType of report
    * \param[in] level is ReportLevel of report
    */
  FileReport(std::string name,
             ReportType type,
             ReportLevel level);
  //@}

  /** @name Destructor */
  //@{
  /**
    * Delete FileReport
    */
  virtual ~FileReport();
  //@}

  /** @name Open methods */
  //@{
  /**
    * Open file
    * \param[in] name is name of file to open
    * \return indicator of success (true) or failure (false)
    */
  virtual bool open(const char* name);
  //@}

 protected:
  /** @name Implementation methods */
  //@{
  /**
    * Print list implementation
    * \param[in] type is ReportType at which to print
    * \param[in] level is ReportLevel at which to print
    * \param[in] format is formatting string
    * \param[in] lst is list of strings to print
    */
  virtual void printListImplementation(ReportType type,
                                       ReportLevel level,
                                       const char* format,
                                       va_list lst);
  /**
    * Flush buffer implementation
    */
  virtual void flushBufferImplementation();
  /**
    * Close implementation
    */
  virtual void closeImplementation();

 private:
  /** @name Default compiler generated methods
    * (Hidden to avoid implicit creation/calling.)
    */
  //@{
  /**
    * Constructor with no arguments
    */
  FileReport();
  /**
    * Copy constructor
    */
  FileReport(const FileReport&);
  /**
    * Overloaded equals operator
    */
  void operator=(const FileReport&);
  //@}

  /**
   * Private members
   */
  //@{
  FILE* file_;
  //@}

};  // end FileReport

/**
  * StreamReport class
  */
class StreamReport : public Report
{

 public:
  /** @name Constructors */
  //@{
  /**
    * Declare StreamReport
    * \param[in] name is name of report
    * \param[in] type is ReportType of report
    * \param[in] level is ReportLevel of report
    */
  StreamReport(std::string name,
               ReportType type,
               ReportLevel level);
  //@}

  /** @name Destructor */
  //@{
  /**
    * Delete StreamReport
    */
  virtual ~StreamReport() {}
  //@}

  /** @name Set methods */
  //@{
  /**
    * Set stream
    * \param[in] os is ostream at which to set stream for report
    */
  void setStream(std::ostream* os);
  //@}

 protected:
  /** @name Implementation methods */
  //@{
  /**
    * Print list implementation
    * \param[in] type is ReportType at which to print
    * \param[in] level is ReportLevel at which to print
    * \param[in] format is formatting string
    * \param[in] lst is list of strings to print
    */
  virtual void printListImplementation(ReportType type,
                                       ReportLevel level,
                                       const char* format,
                                       va_list lst);
  /**
    * Flush buffer implementation
    */
  virtual void flushBufferImplementation();
  /**
    * Close implementation
    */
  virtual void closeImplementation(){};

 private:
  /** @name Default compiler generated methods
    * (Hidden to avoid implicit creation/calling.)
    */
  //@{
  /**
    * Constructor with no arguments
    */
  StreamReport();
  /**
    * Copy constructor
    */
  StreamReport(const StreamReport&);
  /**
    * Overloaded equals operator
    */
  void operator=(const StreamReport&);
  //@}

  /**
   * Private members
   */
  //@{
  std::ostream* os_;
  char buffer_[32768];
  //@}

};  // end StreamReport

}  // namespace NonOpt

#endif /* __NONOPTREPORTER_HPP__ */
