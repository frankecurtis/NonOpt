// Copyright (C) 2022 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cstdio>
#include <cstring>

#include "NonOptReporter.hpp"

namespace NonOpt
{

//////////////
// Reporter //
//////////////

// Constructor
Reporter::Reporter() {}

// Destructor
Reporter::~Reporter()
{

  // Clear elements
  reports_.clear();

} // end destructor

// printf
void Reporter::printf(ReportType type,
                      ReportLevel level,
                      const char* format,
                      ...) const
{

  // Wrap arguments and pass to printList
  va_list lst;
  va_start(lst, format);
  printList(type, level, format, lst);
  va_end(lst);

} // end printf

// Print list
void Reporter::printList(ReportType type,
                         ReportLevel level,
                         const char* format,
                         va_list lst) const
{

  // Print in every report that accepts (type,level)
  for (int i = 0; i < (int)reports_.size(); i++) {
    if (reports_[i]->isAccepted(type, level)) {
      va_list lstcopy;
      va_copy(lstcopy, lst);
      reports_[i]->printList(type, level, format, lstcopy);
      va_end(lstcopy);
    } // end if
  }   // end for

} // end printList

// Add report
bool Reporter::addReport(const std::shared_ptr<Report> report)
{

  // Push
  reports_.push_back(report);

  // Return
  return true;

} // end addReport

// Add file report
bool Reporter::addFileReport(std::string report_name,
                             std::string file_name,
                             ReportType type,
                             ReportLevel level)
{

  // Declare new FileReport
  std::shared_ptr<FileReport> temp(new FileReport(report_name, type, level));

  // Attempt to open file
  if (temp->open(file_name.c_str()) && addReport(temp)) {
    return true;
  }

  // Return
  return false;

} // end addFileReport

// Get report
std::shared_ptr<Report> Reporter::report(std::string name)
{

  // Initialize return value
  std::shared_ptr<Report> returnValue = nullptr;

  // Look for report
  for (int i = 0; i < (int)reports_.size(); i++) {
    std::shared_ptr<Report> temp = reports_[i];
    if (temp->name().compare(name) == 0) {
      returnValue = temp;
      break;
    } // end if
  }   // end for

  // Return
  return returnValue;

} // end report

// Flush buffer
void Reporter::flushBuffer() const
{

  // Flush all buffers
  for (int i = 0; i < (int)reports_.size(); i++) {
    reports_[i]->flushBuffer();
  }

} // end flushBuffer

// Delete reports
void Reporter::deleteReports()
{

  // Delete all reports
  for (int i = 0; i < (int)reports_.size(); i++) {
    reports_[i]->close();
    reports_[i] = nullptr;
  } // end for

  // Resize
  reports_.resize(0);

} // end deleteReports

////////////
// Report //
////////////

// Check accepted
bool Report::isAccepted(ReportType type,
                        ReportLevel level) const
{

  // Check (type,level)
  if (type_ == R_NL) {
    if (type == type_ && level <= level_) {
      return true;
    }
  } // end if
  else {
    if (type == type_ && level == level_) {
      return true;
    }
  } // end else

  // Return
  return false;

} // end isAccepted

////////////////
// FileReport //
////////////////

// Destructor
FileReport::~FileReport()
{

  // Close file
  if (file_ && file_ != stdout && file_ != stderr) {
    fclose(file_);
  }

  // Set pointer to null
  file_ = nullptr;

} // end destructor

// Open file
bool FileReport::open(const char* name)
{

  // If open already, then close it
  if (file_ &&
      file_ != stdout &&
      file_ != stderr) {
    fclose(file_);
  }

  // Set pointer to null
  file_ = nullptr;

  // Set pointer
  if (strcmp("stdout", name) == 0) {
    file_ = stdout;
    return true;
  } // end if
  else if (strcmp("stderr", name) == 0) {
    file_ = stderr;
    return true;
  } // end else if
  else {
    file_ = fopen(name, "w+");
    if (file_) {
      return true;
    }
  } // end else

  // Return default
  return false;

} // end open

// Print list
void FileReport::printList(ReportType type,
                           ReportLevel level,
                           const char* format,
                           va_list lst)
{

  // Print string list
  if (file_) {
    vfprintf(file_, format, lst);
  }

} // end printList

// Flush buffer
void FileReport::flushBuffer()
{

  // Flush buffer
  if (file_) {
    fflush(file_);
  }

} // end flushBuffer

// Close
void FileReport::close()
{

  // Close file
  if (file_ && file_ != stdout && file_ != stderr) {
    fclose(file_);
  }

  // Set pointer to null
  file_ = nullptr;

} // end close

//////////////////
// StreamReport //
//////////////////

// Print list
void StreamReport::printList(ReportType type,
                             ReportLevel level,
                             const char* format,
                             va_list lst)
{

  // Print string list
  if (os_) {
    vsprintf(buffer_, format, lst);
    *os_ << buffer_;
  }

} // end printList

// Flush buffer
void StreamReport::flushBuffer()
{

  // Flush buffer
  if (os_) {
    *os_ << std::flush;
  }

} // end flushBuffer

} // namespace NonOpt
