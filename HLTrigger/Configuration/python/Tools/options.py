# available "type"s and relative global tags
globalTag = {
  'FULL': 'auto:startup',
  'GRun': 'auto:startup',       # use as default
  'data': 'auto:hltonline',
  'HIon': 'auto:starthi',
}


# type used to store a reference to an L1 menu
class ConnectionL1TMenu(object):
  def __init__(self, value):
    self.override = None
    self.connect  = None

    # extract the connection string and configuration name
    if value:
      if ':' in value:
        self.override = 'L1GtTriggerMenu_%s_mc' % value.rsplit(':', 1)[1]
        self.connect  = value.rsplit(':', 1)[0]
      else:
        self.override = 'L1GtTriggerMenu_%s_mc' % value
        self.connect  = None


# type used to store a reference to an HLT configuration
class ConnectionHLTMenu(object):
  def __init__(self, value):
    self.value = value
    self.db    = None
    self.name  = None
    self.run   = None

    # extract the database and configuration name
    if value:
      if ':' in self.value:
        (db, name) = self.value.split(':')
        if db == 'run':
          self.run  = name
        elif db in ('hltdev', 'orcoff'):
          self.db   = db
          self.name = name
        else:
          raise Exception('Unknown ConfDB database "%s", valid values are "hltdev" (default) and "orcoff")' % db)
      else:
        self.db   = 'hltdev'
        self.name = self.value


# options marked with a (*) only apply when creating a whole process configuration
class HLTProcessOptions(object):
  def __init__(self):
    self.menu       = None        #     hlt menu
    self.name       = None        # (*) if set, override the process name
    self.type       = 'GRun'      #     defines global options for 'GRun', 'HIon' or 'online' menus
    self.data       = True        #     run on data (true) or mc (false)
    self.online     = False       # (*) run online (true) or offline (false)
    self.globaltag  = None        # (*) if set, override the GlobalTag
    self.l1         = None        # (*) if set, override the L1 menu
    self.unprescale = False       # (*) if set, unprescale all paths
    self.open       = False       #     if set, cms.ignore all filters, making all paths run on and accept all events
    self.timing     = False       #     if set, instrument the menu for timing measurements
    self.output     = 'all'       # (*) output 'all', 'minimal' or 'none' output modules
    self.fragment   = False       #     prepare a configuration fragment (true) or a whole process (false)
    self.fastsim    = False       #     prepare a configuration fragment suitable for FastSim


  # convert HLT and L1 menus to a dedicated object representation on the fly
  def __setattr__(self, name, value):
    if name is 'menu' and type(value) is not ConnectionHLTMenu:
      # format 'menu' as needed
      object.__setattr__(self, name, ConnectionHLTMenu(value))
    elif name is 'l1' and type(value) is not ConnectionL1TMenu:
      # format '--l1' as needed
      object.__setattr__(self, name, ConnectionL1TMenu(value))
    elif name is 'fastsim' and value:
      # '--fastsim' implies '--fragment' and '--mc'
      object.__setattr__(self, 'fastsim',    True)
      object.__setattr__(self, 'fragment',   True)
      object.__setattr__(self, 'data',       False)
    elif name is 'open' and value:
      # '--open' implies '--unprescale'
      object.__setattr__(self, 'open',       True)
      object.__setattr__(self, 'unprescale', True)
    elif name is 'timing' and value:
      # '--timing' implies '--no-output'
      object.__setattr__(self, 'timing',     True)
      object.__setattr__(self, 'output',     'none')
    else:
      object.__setattr__(self, name, value)
