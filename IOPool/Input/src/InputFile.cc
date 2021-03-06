/*----------------------------------------------------------------------

Holder for an input TFile.
----------------------------------------------------------------------*/
#include "InputFile.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/TimeOfDay.h"

#include <iomanip>

namespace edm {
  InputFile::InputFile(char const* fileName, char const* msg) : file_(), fileName_(fileName), reportToken_(0) {
    logFileAction(msg, fileName);
    file_.reset(TFile::Open(fileName));
    if(!file_) { 
      return;
    }
    if(file_->IsZombie()) { 
      file_.reset();
      return;
    }
    logFileAction("  Successfully opened file ", fileName);
  }

  InputFile::~InputFile() { 
    Close();
  }

  void
  InputFile::inputFileOpened(std::string const& logicalFileName,
                             std::string const& inputType,
                             std::string const& moduleName,
                             std::string const& label,
                             std::string const& fid,
                             std::vector<std::string> const& branchNames) {
    Service<JobReport> reportSvc;
    reportToken_ = reportSvc->inputFileOpened(fileName_,
                                              logicalFileName,
                                              std::string(),
                                              inputType,
                                              moduleName,
                                              label,
                                              fid,
                                              branchNames);
  }

  void
  InputFile::eventReadFromFile(unsigned int run, unsigned int event) const {
    Service<JobReport> reportSvc;
    reportSvc->eventReadFromFile(reportToken_, run, event);
  }

  void
  InputFile::reportInputRunNumber(unsigned int run) const {
    Service<JobReport> reportSvc;
    reportSvc->reportInputRunNumber(run);
  }

  void
  InputFile::reportInputLumiSection(unsigned int run, unsigned int lumi) const {
    Service<JobReport> reportSvc;
    reportSvc->reportInputLumiSection(run, lumi);
  }

  void
  InputFile::Close() {
    if(file_->IsOpen()) {
      file_->Close();
      try {
        logFileAction("  Closed file ", fileName_.c_str());
        Service<JobReport> reportSvc;
        reportSvc->inputFileClosed(reportToken_);
      } catch(...) {
        // If Close() called in a destructor after an exception throw, the services may no longer be active.
        // Therefore, we catch any new exception.
      }
    }
  }

  void
  InputFile::logFileAction(char const* msg, char const* fileName) const {
    LogAbsolute("fileAction") << std::setprecision(0) << TimeOfDay() << msg << fileName;
    FlushMessageLog();
  }
}
