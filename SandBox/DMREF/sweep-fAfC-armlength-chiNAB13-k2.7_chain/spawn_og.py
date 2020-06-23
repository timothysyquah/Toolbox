#!/usr/bin/env python3

import numpy as np
import os
import pdb
import shutil

class Chain:
    ''' Class that stores information about the a given bottlebrush chain
        This will include:
            - backbone length
            - number of A arms (and their starting and ending positions)
            - number of B arms (and their starting and ending positions)
            - the volume fraction of this chain
    '''
    #def __init__(self, Nbb, Nsc, naarm, agraftstart, agraftend, nbarm, bgraftstart, bgraftend):
    #    self.Nbb = Nbb
    #    self.Nsc = Nsc
    #    self.naarm = naarm
    #    self.agraftstart = agraftstart
    #    self.agraftend = agraftend
    #    self.nbarm = nbarm
    #    self.bgraftstart = bgraftstart
    #    self.bgraftend = bgraftend

    def __init__(self,fAfC, delta, narmtotal, Nsc):
        
        
        fA = fAfC
        fC = fAfC
        fB = 1.0 - fA - fC

        self.nblock = 3 
        self.nameblock = ['A','B','C']
        self.fblock = [ fA, fB, fC ]
        assert ( np.isclose(sum(self.fblock),1)), f"block fractions must sum to 1 {self.fblock} {sum(self.fblock)}"

        self.narmblock = [None]*self.nblock
        self.Nscblock = [None]*self.nblock
        for i in range(self.nblock):
            narm = narmtotal*self.fblock[i]
            assert(np.isclose(narm,np.round(narm))), f"narm of {i}th block must be an integer {narm}"
            narm = np.round(narm)
            self.narmblock[i] = narm
            self.Nscblock[i] = Nsc

        self.Nbb=delta*(narmtotal-1)

        self.arm_position_type = 'start-end'

        # determine graft start and end positions for each block
        self.graftstart = [None]*self.nblock
        self.graftend = [None]*self.nblock
        for i in range(self.nblock):
            if i == 0:
                self.graftstart[i] =0.0;
            else:   
                self.graftstart[i] = self.graftend[i-1]+delta

            if i == (self.nblock - 1):
                self.graftend[i] = self.Nbb
            else:
                self.graftend[i] = self.graftstart[i] + (self.narmblock[i]-1) * delta


    def write_to_file(self, filehandle, ichain):
        f = filehandle
        
        # check if all arms are zero length -- if so its just a linear chain
        islinearchain=True
        for b in range(self.nblock): 
          if self.Nscblock[b] != 0:
            islinearchain=False


        if islinearchain:
          f.write(f'''
            chain{ichain} {{
                label  = ABC
                architecture = linear
                statistics = gaussian
                length       = {self.Nbb : 0.5f}
                nblocks      = {self.nblock}
                blockspecies = 1 2 3 
                blockfractions = ''')
          for i in range(self.nblock):
              f.write(f'{self.fblock[i] : 0.5f} ')
          f.write('''
                  }\n''')

        else: # bottlebrush
          narmtypes = self.nblock
          #  write chain comb
          f.write(f'''
            chain{ichain} {{
                label  = bottlebrush{ichain}
                architecture = comb
                NumSideArmTypes = {narmtypes}
            ''')

          f.write(f'''
                backbone {{ 
                  statistics   = gaussian
                  length       = {self.Nbb : 0.5f}
                  nblocks      = {self.nblock}
                  blockspecies = 1 2 3 
                  blockfractions = ''')
          for i in range(self.nblock):
              f.write(f'{self.fblock[i] : 0.5f} ')
          f.write('''
                  }\n''')

          for i in range(self.nblock):
              Nsc = self.Nscblock[i]
              narm = self.narmblock[i]
              assert(self.arm_position_type == 'start-end')

              graftstart = self.graftstart[i]
              graftend = self.graftend[i]
              blocktype = i + 1

              f.write(f'''
                sidearmtype{i+1} {{
                  statistics   = gaussian
                  length       = {Nsc : 0.5f}''')

              # write block info
              f.write(f'''
                  nblocks      = 1
                  blockspecies = {blocktype}
                  NumArms               =  {narm}
                  BackboneGraftingStart = {graftstart : 0.5f}
                  BackboneGraftingEnd   = {graftend : 0.5f}
                  }}''')
          f.write('''
            }''')


if __name__ == "__main__":
    #-----------------------------------------------------
    # Parameters 
   
    #PhaseList      = ['DIS', 'LAM', 'ACYL', 'NaCl', 'CsCl', 'AGYR32','ADIA32', 'CsCl64', 'AGYR64']
    #InitialL0Guess = [3.0, 3.0, 3.5, 4.0, 3.5, 5.0,5.0, 3.5,5.0]
    #NThreads       = [  1,   1,   1,   1, 1, -1,-1, -1,-1]

    #PhaseList      = ['DIS', 'CsCl','NaCl'] 
    #InitialL0Guess = [3.0, 3.5,3.0]
    #NThreads       = [ 1, 1,1]

    #PhaseList      = ['AGYR32', 'ADIA32', 'ACYL','DIS'] 
    #InitialL0Guess = [5.0, 5.0,3.5,3.0]
    #NThreads       = [ -1, -1,1,1]

    PhaseList      = ['AGYR32', 'ADIA32', 'ACYL'] 
    SpaceGroup     = ['I4_1.32','Fd-3m:1','p4mm']
    InitialL0Guess = [5.0, 5.0,3.5]
    NThreads       = [ -1, -1,1]
    
    #PhaseList      = ['AGYRC'] 
    #InitialL0Guess = [6.0]
    #NThreads       = [-1]

    #PhaseList      = ['ADIAC64'] 
    #InitialL0Guess = [6.0]
    #NThreads       = [-1]

    Nscstart = 0.040
    Nscend = 0.040
    Nscdelta = 0.020
    fAfCstart = 0.10 # 0.2
    fAfCend = 0.30 # 0.2
    fAfCdelta = 0.02 # 0.01

    #=======================================================
    # FIXED Parameters - probably shouldn't be changing
    delta=0.01 #~6k  #3k - 0.025 #0.01 # sidechain spacing
    DS=0.002
    narmtotal=100
    chiNAB = 13.0 / (delta*(narmtotal-1)) # note: normalizing by backbone length since in linear chain this is 1
    k = 2.7
    chiNBC = chiNAB
    chiNAC = chiNAB*k
    #=======================================================

    # populate lists
    if Nscend == Nscstart:
      nNsc = 1
    else:
      nNsc = (Nscend - Nscstart) / Nscdelta 
      assert (np.isclose(nNsc,round(nNsc,3))), f"nNsc {nNsc} is not an int. Change bounds!"
      nNsc = int(nNsc) + 1 # include both end points
    Nscs = np.linspace(Nscstart,Nscend,nNsc)
    if fAfCend == fAfCstart:
      nfAfC = 1
    else:
      nfAfC = (fAfCend - fAfCstart) / fAfCdelta 
      assert (np.isclose(nfAfC,round(nfAfC,3))), f"nfAfC {nfAfC} is not an int. Change bounds!"
      nfAfC = round(nfAfC) + 1 # include both end points
    fAfCs = np.linspace(fAfCstart,fAfCend,nfAfC)
    print(Nscs)
    print(fAfCs)

    try: 
        assert(len(PhaseList) == len(InitialL0Guess) == len(NThreads))
    except:
        if type(NThreads) == int and type(InitialL0Guess) == float and type(PhaseList) == str:
            pass
        else:
            raise RuntimeError("Lengths of PhasesList, InitialL0Guess and NThreads not equal! (%d != %d = %d)" % (len(PhaseList),len(InitialL0Guess),len(NThreads)) )

    IDIR = os.getcwd()

    for fAfC in fAfCs:
        for Nsc in Nscs:
            chain = Chain(fAfC,delta,narmtotal,Nsc)

            for iphase,phase in enumerate(PhaseList):
                 
                WDIR=f'Nsc{Nsc:0.3f}/fAfC{fAfC:0.2f}/{phase}Phase/'

                print (WDIR)
                if os.path.exists(WDIR):
                    print("{} exists...skipping.".format(WDIR)) 
                    continue
                else:
                    os.makedirs(WDIR)

            
                # now write *.in file
                infile=f'{phase}.in'
                outfile=f'{phase}.out'
                with open('%s/%s' % (WDIR,infile), 'w') as fout:
                    with open('SEEDS/%s' % infile,'r') as f:
                        for line in f:
                            if 'chains {' in line:
                                line += "\tnchains = 1\n"
                                line += "\tcontourds = %f\n" % DS
                                line += "\tdiffusermethod = SOS\n"

                            #line = line.replace('__NPW__', f'{npw}')
                            #line = line.replace('__DT__',  f'{DT : 0.5f}')
                            line = line.replace('__chiN12__',f'{chiNAB: 0.5f}')
                            line = line.replace('__chiN13__',f'{chiNAC: 0.5f}')
                            line = line.replace('__chiN23__',f'{chiNBC: 0.5f}')
                            line = line.replace('__L0__',  f'{InitialL0Guess[iphase]: 0.5f}')

                            fout.write(line)

                            if 'chains {' in line:
                                chain.write_to_file(fout, 1)

                # now edit submit file
                if NThreads[iphase] == 1:
                    submitfile='SEEDS/submit.sh'
                elif NThreads[iphase] == -1:
                    if phase != 'SIGMA':
                        submitfile='SEEDS/submitGPU.sh'
                    else:
                        submitfile='SEEDS/submitGPU-P100.sh'
                else:
                    submitfile='SEEDS/submitPLL.sh'

                with open('%s/submit.sh' % WDIR, 'w') as fout:
                    with open(submitfile,'r') as f:
                        for line in f:
                            #line = line.replace('__JOBNAME__',WDIR)
                            line = line.replace('__INFILE__',f'{phase}.in')
                            line = line.replace('__OUTFILE__',f'{phase}.out')
                            if NThreads[iphase] > 1:
                                line = line.replace('__NTHREADS__',WDIR)
                            fout.write(line)

                
                # copy seed file if it exists
                seedfieldsfile=f'SEEDS/{phase}_fields.in'
                if os.path.exists(seedfieldsfile):
                   shutil.copyfile(seedfieldsfile,'%s/fields.in' % WDIR)
                    

                os.chdir(WDIR)

#                import subprocess
#                cmd="qsub submit.sh"
                #cmd="sbatch submit.sh"
#                subprocess.call(cmd.split())

                os.chdir(IDIR)


