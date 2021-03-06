#ifndef FUEVENTPROCESSOR_H
#define FUEVENTPROCESSOR_H 1

#include "EventFilter/Utilities/interface/StateMachine.h"
#include "EventFilter/Utilities/interface/RunBase.h"
#include "EventFilter/Utilities/interface/Css.h"
#include "EventFilter/Utilities/interface/Exception.h"
#include "EventFilter/Utilities/interface/SquidNet.h"
#include "EventFilter/Utilities/interface/Vulture.h"

#include "EventFilter/Utilities/interface/MasterQueue.h"
#include "EventFilter/Utilities/interface/SlaveQueue.h"
#include "SubProcess.h"
#include "FWEPWrapper.h"

#include "DataFormats/Provenance/interface/ModuleDescription.h"

#include "xdaq/Application.h"
#include "xdaq/NamespaceURI.h"

#include "xdata/String.h"
#include "xdata/Integer.h"
#include "xdata/Boolean.h"
#include "xdata/Vector.h"
#include "xdata/UnsignedInteger32.h"
#include "xdata/ActionListener.h"
#include "xdata/InfoSpaceFactory.h"

#include "xgi/Input.h"
#include "xgi/Output.h"
#include "xgi/exception/Exception.h"

#include <sys/time.h>
#include <pthread.h>

#include <list>
#include <vector>
#include <map>

namespace evf
{

  /* to be filled in with summary from paths */
  struct filter {

  };
  class CPUStat;
  class RateStat;
  namespace internal{
    
    class MyCgi : public xgi::Input{
    public:
      MyCgi() : xgi::Input("",0){}
      //      MyCgi(xgi::Input &b) : xgi::Input("",0) {environment_ = b.environment_;}
      std::map<std::string, std::string, std::less<std::string> > &getEnvironment(){return environment_;}
    };
  }
  class FUEventProcessor : public xdaq::Application,
			   public xdata::ActionListener
  {
  public:
    //
    // construction/destruction
    //
    XDAQ_INSTANTIATOR();
    FUEventProcessor(xdaq::ApplicationStub *s);
    virtual ~FUEventProcessor();
    

    //
    // member functions
    //

    // work loop functions to be executed during transitional states (async)
    bool configuring(toolbox::task::WorkLoop* wl);
    bool enabling(toolbox::task::WorkLoop* wl);
    bool stopping(toolbox::task::WorkLoop* wl);
    bool halting(toolbox::task::WorkLoop* wl);

    // fsm soap command callback
    xoap::MessageReference fsmCallback(xoap::MessageReference msg)
      throw (xoap::exception::Exception);
        
    // xdata:ActionListener interface
    void actionPerformed(xdata::Event& e);
    
    // trigger report related helper functions
    //    std::string triggerReportToString(const edm::TriggerReport& tr);
    //    void triggerReportToTable(const edm::TriggerReport& tr);
    //    void        printTriggerReport(const edm::TriggerReport& tr);

    // HyperDAQ related functions
    void defaultWebPage(xgi::Input *in,xgi::Output *out)
      throw(xgi::exception::Exception);
    void spotlightWebPage(xgi::Input *,xgi::Output *)
      throw(xgi::exception::Exception);
    void scalersWeb(xgi::Input *,xgi::Output *)
      throw(xgi::exception::Exception);
    void pathNames(xgi::Input *,xgi::Output *)
      throw(xgi::exception::Exception);
    void css(xgi::Input *in,xgi::Output *out) throw (xgi::exception::Exception)
    {
      css_.css(in,out);
    }

    void subWeb(xgi::Input *in,xgi::Output *out);
    void moduleWeb(xgi::Input *in,xgi::Output *out){evtProcessor_.moduleWeb(in,out);}
    void serviceWeb(xgi::Input *in,xgi::Output *out){evtProcessor_.serviceWeb(in,out);}
    void microState(xgi::Input *in,xgi::Output *out);
    void updater(xgi::Input *in,xgi::Output *out);
    void procStat(xgi::Input *in,xgi::Output *out);
    void sendMessageOverMonitorQueue(MsgBuf &);

  private:

    void attachDqmToShm()   throw (evf::Exception);
    void detachDqmFromShm() throw (evf::Exception);

    std::string logsAsString();
    void localLog(std::string);

    // MPEP functions
    void startSupervisorLoop();
    void startReceivingLoop();
    void startReceivingMonitorLoop();
    // calculate scalers information in separate thread
    void startScalersWorkLoop() throw (evf::Exception);
    bool scalers(toolbox::task::WorkLoop* wl);
    void startSummarizeWorkLoop() throw (evf::Exception);
    bool summarize(toolbox::task::WorkLoop* wl);

    bool receiving(toolbox::task::WorkLoop* wl);
    bool receivingAndMonitor(toolbox::task::WorkLoop* wl);
    bool supervisor(toolbox::task::WorkLoop* wl);
    bool enableCommon();
    bool enableClassic();
    //    void enableMPEPMaster();
    bool enableMPEPSlave();
    bool stopClassic();
    void stopSlavesAndAcknowledge();
    void makeStaticInfo();

    //
    // member data
    //
    
    // finite state machine
    evf::StateMachine                fsm_;
    
    // logger
    Logger                           log_;

    // edm event processor
    FWEPWrapper                      evtProcessor_;
    
    // parameters published to XDAQ info space(s)
    xdata::String                    url_;
    xdata::String                    class_;
    xdata::UnsignedInteger32         instance_;
    xdata::UnsignedInteger32         runNumber_;
    xdata::Boolean                   epInitialized_; 
    xdata::String                    configString_;
    std::string                      configuration_;

    xdata::Boolean                   outPut_;

    xdata::Boolean                   autoRestartSlaves_;
    xdata::UnsignedInteger32         slaveRestartDelaySecs_;

    xdata::Boolean                   hasShMem_;
    xdata::Boolean                   hasPrescaleService_;
    xdata::Boolean                   hasModuleWebRegistry_;
    xdata::Boolean                   hasServiceWebRegistry_;
    xdata::Boolean                   isRunNumberSetter_;
    xdata::Boolean                   iDieStatisticsGathering_;
    bool                             outprev_;
    
    // application identifier
    std::string                      sourceId_;

    // flashlist variables, squids
    xdata::Boolean                   squidPresent_; 

    // behavior on error - configurable
    
    xdata::Boolean                   exitOnError_; 

    // HyperDAQ related
    Css                              css_;

    // Misc
    std::string                      reasonForFailedState_;

    SquidNet                         squidnet_;
    std::vector<std::string>         logRing_;
    unsigned int                     logRingIndex_;
    static const unsigned int        logRingSize_ = 50;
    bool                             logWrap_;

    xdata::UnsignedInteger32         nbSubProcesses_;
    xdata::UnsignedInteger32         nbSubProcessesReporting_;
    std::vector<SubProcess>          subs_;
    unsigned int                     nblive_; 
    unsigned int                     nbdead_; 

    unsigned int                     nbTotalDQM_;

    // workloop / action signature for message passing
    toolbox::task::WorkLoop         *wlReceiving_;      
    toolbox::task::ActionSignature  *asReceiveMsgAndExecute_;
    bool                             receiving_;
    toolbox::task::WorkLoop         *wlReceivingMonitor_;      
    toolbox::task::ActionSignature  *asReceiveMsgAndRead_;
    bool                             receivingM_;
    SubProcess*                      myProcess_;
    toolbox::task::WorkLoop         *wlSupervising_;      
    toolbox::task::ActionSignature  *asSupervisor_;
    bool                             supervising_;

    xdata::InfoSpace*                monitorInfoSpace_;
    xdata::InfoSpace*                monitorLegendaInfoSpace_;
    xdata::InfoSpace*                applicationInfoSpace_;
    pthread_mutex_t                  stop_lock_;
    pthread_mutex_t                  start_lock_;
    pthread_mutex_t                  pickup_lock_;
    std::string                      updaterStatic_;
    xdata::Serializable             *nbProcessed;
    xdata::Serializable             *nbAccepted;

    // flahslist variables, scalers
    xdata::InfoSpace                *scalersInfoSpace_;
    xdata::InfoSpace                *scalersLegendaInfoSpace_;

    //scalers workloop
    toolbox::task::WorkLoop         *wlScalers_;      
    toolbox::task::ActionSignature  *asScalers_;
    bool                             wlScalersActive_;
    unsigned int                     scalersUpdates_;

    //summarize workloop
    toolbox::task::WorkLoop         *wlSummarize_;      
    toolbox::task::ActionSignature  *asSummarize_;
    bool                             wlSummarizeActive_;
    int                              anonymousPipe_[2];
    xdata::Vector<xdata::Integer>    spMStates_;
    xdata::Vector<xdata::Integer>    spmStates_;
    xdata::UnsignedInteger32         superSleepSec_; 
    std::list<std::string>           names_;
    xdata::String                    iDieUrl_;
    Vulture                         *vulture_;
    pid_t			     vp_;
    CPUStat                         *cpustat_;
    RateStat                        *ratestat_;
  };
  
} // namespace evf


#endif
