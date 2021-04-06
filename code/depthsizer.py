#!/usr/bin/python

# See below for name and description
# Copyright (C) 2016 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
#  
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not, write to 
# the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Author contact: <seqsuite@gmail.com> / School of Biotechnology and Biomolecular Sciences, UNSW, Sydney, Australia.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       depthsizer
Description:  Read-depth based genome size prediction
Version:      1.0.0
Last Edit:    06/04/21
Copyright (C) 2021  Richard J. Edwards - See source code for GNU License Notice

Function:
    DepthSizer is a wrapper for the genome size estimate methods of Diploidocus. DepthSizer needs a genome assembly
    (fasta format, `seqin=FILE`), a set of long read (ONT, PacBio or HiFi) data for the assembly (`reads=FILELIST` and
    `readtype=LIST`) (or `readbp=INT`), and a BUSCO full table of results (`busco=TSVFILE`).

    DepthSizer works on the principle that `Complete` BUSCO genes should represent predominantly single copy (diploid
    read depth) regions along with some ooor quality and/or repeat regions. Assembly artefacts and collapsed repeats etc.
    are predicted to deviate from diploid read depth in an inconsistent manner. Therefore, even if less than half the
    region is actually diploid coverage, the **modal** read depth is expected to represent the actual single copy
    read depth.

    DepthSizer uses `samtools mpileup` (or `samtools depth` if `quickdepth=T`) to calculate the per-base read depth
    and extracts the modal read depth for each BUSCO gene along with the overall modal read depth for all gene
    regions. Genome size is then estimated based on a crude calculation using the total combined sequencing length.
    This will be caculated from `reads=FILELIST` unless provided with `readbp=INT`.

    BUSCO single-copy genes are parsed from a BUSCO full results table, given by `busco=TSVFILE` (default
    `full_table_$BASEFILE.busco.tsv`). This can be replaced with any table with the fields:
    ['BuscoID','Status','Contig','Start','End','Score','Length']. Entries are reduced to those with `Status` = `Complete`
    and the `Contig`, `Start` and `End` fields are used to define the regions that should be predominantly single copy.
    Output from BUSCOMP is also compatible with DepthSizer.

    **NOTE:** The current genome size prediction appears to be an over-estimate. There is currently no adjustment for
    contamination. The `mapadjust` option attemtps to correct for read mapping and imbalanced insertion:deletion ratios
    etc. but has not been extensively tested.

    ---

Dependencies:

    Unless `bam=FILE` is given, [minimap2](https://github.com/lh3/minimap2) must be installed and either added to the
    environment `$PATH` or given to DepthSizer with the `minimap2=PROG` setting, and [samtools](http://www.htslib.org/)
    needs to be installed. Unless `depdensity=F`, R will also need be installed.

    To generate documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

        export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

    For DepthSizer documentation, run with `dochtml=T` and read the `*.docs.html` file generated.

Commandline:
    ### ~ Main DepthSizer run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence assembly [None]
    basefile=FILE   : Root of output file names [gapspanner or $SEQIN basefile]
    summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
    genomesize=INT  : Haploid genome size (bp) [0]
    bam=FILE        : BAM file of long reads mapped onto assembly [$BASEFILE.bam]
    reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    dochtml=T/F     : Generate HTML DepthSizer documentation (*.docs.html) instead of main run [False]
    tmpdir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]
    ### ~ Genome size prediction options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    busco=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
    readbp=INT      : Total combined read length for depth calculations (over-rides reads=FILELIST) []
    quickdepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]
    depdensity=T/F  : Whether to use the BUSCO depth density profile in place of modal depth [True]
    mapadjust=T/F   : Whether to adjust predicted genome size based on read length:mapping ratio [True]
    ### ~ Forking options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    forks=X         : Number of parallel sequences to process at once [0]
    killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
    killmain=T/F    : Whether to kill main thread rather than individual forks when killforks reached. [False]
    logfork=T/F     : Whether to log forking in main log [False]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import diploidocus
import rje, rje_obj, rje_rmd
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 1.0.0 - Initial working version of DepthSizer based on Diploidocus v0.16.2.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [N] : Add REST outputs to restSetup() and restOutputOrder()
    # [Y] : Add to SLiMSuite or SeqSuite.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('DepthSizer', '1.0.0', 'April 2021', '2021')
    description = 'Read-depth based genome size prediction'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_obj.zen()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cmd_help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if cmd_help > 0:
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Major Problem with cmdHelp()')
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program when called from commandline.
    '''
    Basic Setup of Program when called from commandline:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:### ~ [1] ~ Initial Command Setup & Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        info = makeInfo()                                   # Sets up Info object with program details
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('%s v%s' % (info.program,info.version)); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['description','-description','--description']: rje.printf('%s: %s' % (info.program,info.description)); sys.exit(0)
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)   # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2 
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Problem during initial setup.'); raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: DepthSizer Class                                                                                       #
#########################################################################################################################
class DepthSizer(rje_obj.RJE_Object):
    '''
    DepthSizer Class. Author: Rich Edwards (2021).

    Str:str
    
    Bool:boolean

    Int:integer

    Num:float

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = []
        self.boollist = ['DocHTML']
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({})
        self.setInt({})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._forkCmd(cmd)  # Delete if no forking
                ### Class Options (No need for arg if arg = att.lower()) ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['DocHTML'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''
        # DepthSizer: Read-depth based genome size prediction

        DepthSizer is a wrapper for the genome size estimate methods of Diploidocus. DepthSizer needs a genome assembly
        (fasta format, `seqin=FILE`), a set of long read (ONT, PacBio or HiFi) data for the assembly (`reads=FILELIST` and
        `readtype=LIST`) (or `readbp=INT`), and a BUSCO full table of results (`busco=TSVFILE`).

        DepthSizer works on the principle that `Complete` BUSCO genes should represent predominantly single copy (diploid
        read depth) regions along with some ooor quality and/or repeat regions. Assembly artefacts and collapsed repeats etc.
        are predicted to deviate from diploid read depth in an inconsistent manner. Therefore, even if less than half the
        region is actually diploid coverage, the **modal** read depth is expected to represent the actual single copy
        read depth.

        DepthSizer uses `samtools mpileup` (or `samtools depth` if `quickdepth=T`) to calculate the per-base read depth
        and extracts the modal read depth for each BUSCO gene along with the overall modal read depth for all gene
        regions. Genome size is then estimated based on a crude calculation using the total combined sequencing length.
        This will be caculated from `reads=FILELIST` unless provided with `readbp=INT`.

        BUSCO single-copy genes are parsed from a BUSCO full results table, given by `busco=TSVFILE` (default
        `full_table_$BASEFILE.busco.tsv`). This can be replaced with any table with the fields:
        ['BuscoID','Status','Contig','Start','End','Score','Length']. Entries are reduced to those with `Status` = `Complete`
        and the `Contig`, `Start` and `End` fields are used to define the regions that should be predominantly single copy.
        Output from BUSCOMP is also compatible with DepthSizer.

        **NOTE:** The current genome size prediction appears to be an over-estimate. There is currently no adjustment for
        contamination. The `mapadjust` option attemtps to correct for read mapping and imbalanced insertion:deletion ratios
        etc. but has not been extensively tested.

        ---

        # Running DepthSizer

        DepthSizer is written in Python 2.x and can be run directly from the commandline:

            python $CODEPATH/depthsizer.py [OPTIONS]

        If running as part of [SLiMSuite](http://slimsuite.blogspot.com/), `$CODEPATH` will be the SLiMSuite `tools/`
        directory. If running from the standalone [DepthSizer git repo](https://github.com/slimsuite/depthsizer), `$CODEPATH`
        will be the path the to `code/` directory. Please see details in the [DepthSizer git repo](https://github.com/slimsuite/depthsizer)
        for running on example data.

        ## Dependencies

        Unless `bam=FILE` is given, [minimap2](https://github.com/lh3/minimap2) must be installed and either added to the
        environment `$PATH` or given to DepthSizer with the `minimap2=PROG` setting, and [samtools](http://www.htslib.org/)
        needs to be installed. Unless `depdensity=F`, R will also need be installed.

        To generate documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

            export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

        For DepthSizer documentation, run with `dochtml=T` and read the `*.docs.html` file generated.

        ## Commandline options

        A list of commandline options can be generated at run-time using the `-h` or `help` flags. Please see the general
        [SLiMSuite documentation](http://slimsuite.blogspot.com/2013/08/command-line-options.html) for details of how to
        use commandline options, including setting default values with **INI files**.

        ```
        ### ~ Main DepthSizer run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        seqin=FILE      : Input sequence assembly [None]
        basefile=FILE   : Root of output file names [depthsizer or $SEQIN basefile]
        summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
        genomesize=INT  : Haploid genome size (bp) [0]
        bam=FILE        : BAM file of long reads mapped onto assembly [$BASEFILE.bam]
        reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
        readtype=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
        dochtml=T/F     : Generate HTML DepthSizer documentation (*.docs.html) instead of main run [False]
        tmpdir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]
        ### ~ Genome size prediction options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        busco=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
        readbp=INT      : Total combined read length for depth calculations (over-rides reads=FILELIST) []
        quickdepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]
        depdensity=T/F  : Whether to use the BUSCO depth density profile in place of modal depth [True]
        mapadjust=T/F   : Whether to adjust predicted genome size based on read length:mapping ratio [True]
        ### ~ Forking options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        forks=X         : Number of parallel sequences to process at once [0]
        killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
        killmain=T/F    : Whether to kill main thread rather than individual forks when killforks reached. [False]
        logfork=T/F     : Whether to log forking in main log [False]
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        # DepthSizer workflow and options

        The main inputs for DepthSizer genome size prediction are:

        * `seqin=FILE` : Input sequence assembly to tidy [Required].
        * `reads=FILELIST`: List of fasta/fastq files containing long reads. Wildcard allowed. Can be gzipped. For a single run (not cycling), a BAM file can be supplied instead with `bam=FILE`. (This will be preferentially used if found, and defaults to `$BASEFILE.bam`.) Read types (pb/ont/hifi) for each file are set with `readtype=LIST`, which will be cycled if shorter (default=`ont`). Optionally, the pre-calculated total read length can be provided with `readbp=INT` and/or the pre-calculated (haploid) genome size can be provided with `genomesize=INT`.
        * `readtype=LIST` : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
        * `busco=TSVFILE` : BUSCO full table [`full_table_$BASEFILE.busco.tsv`] used for calculating single copy ("diploid") read depth.

        ## Step 1: BAM file (read mapping)

        The first step is to generate a BAM file by mapping `reads` on to `seqin` using [minimap2](https://github.com/lh3/minimap2).
        A pre-generated BAM file can be given instead using `bam=FILE`. There should be no secondary mapping of reads, as
        these will inflate read depths, so filter these out if they were allowed during mapping. Similarly, the BAM file
        should not contain unmapped reads.

        ## Step 2: BUSCO(MP) results

        DepthSizer works on the principle that `Complete` BUSCO genes should represent predominantly single copy (diploid read
        depth) regions along with some poor-quality and/or repeat regions. Assembly artefacts and collapsed repeats etc.
        are predicted to deviate from diploid read depth in an inconsistent manner. Therefore, even if less than half the
        region is actually diploid coverage, the **modal** read depth is expected to represent the actual single copy
        read depth.

        BUSCO single-copy genes are parsed from a BUSCO full results table, given by `busco=TSVFILE` (default
        `full_table_$BASEFILE.busco.tsv`). This can be replaced with any table with the fields:
        ['BuscoID','Status','Contig','Start','End','Score','Length']. Entries are reduced to those with `Status` = `Complete`
        and the `Contig`, `Start` and `End` fields are used to define the regions that should be predominantly single copy.
        [BUSCOMP](https://github.com/slimsuite/buscomp) v0.10.0 and above will generate a `*.complete.tsv` file that can
        be used in place of BUSCO results. This can enable rapid re-annotation of BUSCO genes following, for example,
        vector trimming with [Diploidocus](https://github.com/slimsuite/diploidocus).

        ## Step 3: Total read volume

        Genome size is estimated based on a crude calculation using the total combined sequencing length. This will be
        calculated from `reads=FILELIST` unless provided with `readbp=INT`.

        The `mapadjust` option attempts to correct for read mapping and imbalanced insertion:deletion ratios
        etc. but has not been extensively tested. This uses `samtools coverage` to estimate the total number of bases
        mapped onto the assembly (assembly bases with coverage x average depth) and `samtools fasta` to extract all the
        mapped reads. The ratio
        of the mapped read bases to the summed length of mapped reads is then calculated and used to adjust the total
        combined read length. This is an attempt to adjust for losses in mapped reads at sequence ends and/or gaps (for
        fragmented assemblies) and any sequencing error bias towards insertions or deletions. This adjustment can be
        switched off with `mapadjust=F`.

        **NOTE:** when running through Diploidocus, `mapadjust` is set to `False` by default for consistency with older runs.

        ## Step 4: Single-copy read depth

        DepthSizer uses `samtools mpileup` (or `samtools depth` if `quickdepth=T`) to calculate the per-base read depth
        and extracts the modal read depth for each BUSCO gene along with the overall modal read depth for all BUSCO gene
        regions. Unless `depdensity=F`, the R function density() is used to find the peak of the overall BUSCO read depth
        and this is used in place of the modal read depth. DepthSizer will report predicted genome sizes using all three
        methods (see below). Preliminary analyses indicate that the density-derived peak depth is most accurate, but this
        may change with further testing. Each depth estimate will be output as `#MODE` entries in the log.

        BUSCO genes are used to self-validate the single-copy read depth estimate and the estimated copy number of each
        gene is output to `*.buscoX.tdt`, where `X` is the depth method used (`mpileup` or `depth`) with mean values
        output as `#CNV` log entries.

        ## Step 5: Genome size prediction

        The final genome size is predicted based on the total (adjusted) combined sequencing length and the single-copy
        read depth, as: `readbp`/`scdepth`.

        Size predictions will be output for all single-copy read depth estimates calculated. These are found as `#GSIZE`
        entries in the log file.

        **NOTE:** There is currently no adjustment for non-nuclear DNA or contamination. The current unadjusted genome size
        prediction appears to be an over-estimate. The `mapadjust` setting is still under development and may be prone to
        issues for highly fragmented/gappy genome assemblies. Extreme mapadjust ratios should be treated with caution and
        may indicate problems with the assembly and/or source data.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DocHTML'): return rje_rmd.docHTML(self)
            dipobj = diploidocus.Diploidocus(self.log,['dna=T','summarise=T','mapadjust=T']+self.cmd_list+['runmode=gensize','diploidocus=F'])
            if not dipobj.getStrLC('SeqIn'):
                raise ValueError('seqin=FILE must be set')
            if not rje.exists(dipobj.getStr('SeqIn')):
                self.errorLog('seqin=FILE "{0}" not found!',format(dipobj.getStr('SeqIn')),printerror=False)
                raise IOError()
            #if dipobj.baseFile() == 'diploidocus': dipobj.baseFile('gapsizer')
            if not dipobj.setup(): return False
            self.baseFile(dipobj.baseFile())
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dipobj.seqinObj(summarise=False)
            dipobj.genomeSize(makebam=True)
            dipobj.db('gensize').saveToFile(append=self.getBool('Append'))
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
### End of SECTION II: DepthSizer Class                                                                                 #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################

#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return  
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: DepthSizer(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
