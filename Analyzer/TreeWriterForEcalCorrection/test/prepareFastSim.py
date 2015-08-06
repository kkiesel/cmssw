#!/bin/env python
import os
import re
import glob
import string

from optparse import OptionParser


def createFromTemplate( d, templateFilename, outFilename  ):
   o = open(outFilename,"w+")
   data = open(templateFilename).read()
   tem=string.Template(data)
   data=tem.safe_substitute(d)
   o.write(data)
   o.close()

def prepareFastSim( isFull ):

    fname = "full" if isFull else "fast"

    cfgTemplateName = "SingleElectron_{}_template_cfi_GEN_SIM.py".format( fname )
    crabTemplateName = "crab.cfg"


    tmpFolder = "crab"
    os.mkdir( tmpFolder )
    os.chdir( tmpFolder )

    for energy in range( 10, 100, 10 ) + [ 100, 150, 200, 300, 400, 600, 800, 1000 ]:

        baseName = "SingleElectronGun_E{}_{}"
        name = baseName.format( energy, fname )

        crabReplacementDict = {
            "PSET": cfgTemplateName,
            "NAME": name
        }
        createFromTemplate( crabReplacementDict, "../"+crabTemplateName, crabTemplateName )

        cfgReplacementDict = {
            "NAME": name,
            "ENERGY": energy
        }
        createFromTemplate( cfgReplacementDict, "../"+cfgTemplateName, cfgTemplateName )

        os.system( "crab -create -submit" )

if (__name__ == "__main__"):
    # use option parser to allow verbose mode
    parser = OptionParser()
    parser.add_option( "--full", action="store_true", default=False )
    (opts, args) = parser.parse_args()

    prepareFastSim( opts.full )
