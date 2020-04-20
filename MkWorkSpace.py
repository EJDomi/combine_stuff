#!/usr/bin/env python

from __future__ import print_function
import os, sys, re, imp, ROOT, copy
import time
from collections import OrderedDict
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit.so')

debug = False

datacard_template='''----------------------------------------------------------------------------------------------------------------------------------
imax 1 number of bins
jmax 5 number of processes minus 1
kmax * number of nuisance parameters
----------------------------------------------------------------------------------------------------------------------------------
shapes data_obs          REG_NAME  WS w:data_obs_MPE_REG
shapes PROC_MASS         REG_NAME  WS w:PROC_MASS_MPE_REG w:PROC_MASS_MPE_REG_$SYSTEMATIC
shapes TTJets2017        REG_NAME  WS w:TTJets2017_MPE_REG w:TTJets2017_MPE_REG_$SYSTEMATIC
shapes WJets2017         REG_NAME  WS w:WJets2017_MPE_REG w:WJets2017_MPE_REG_$SYSTEMATIC
shapes ZJets2017         REG_NAME  WS w:ZJets2017_MPE_REG w:ZJets_2017_MPE_REG_$SYSTEMATIC
shapes VV2017            REG_NAME  WS w:VV2017_MPE_REG w:VV2017_MPE_REG_$SYSTEMATIC
shapes DY2017            REG_NAME  WS w:DY2017_MPE_REG w:DY2017_MPE_REG_$SYSTEMATIC
----------------------------------------------------------------------------------------------------------------------------------
bin                REG_NAME      
observation        -1.0    
----------------------------------------------------------------------------------------------------------------------------------
bin                REG_NAME        REG_NAME      REG_NAME      REG_NAME     REG_NAME      REG_NAME
process            PROC_MASS      TTJets2017     WJets2017    ZJets2017      VV2017        DY2017
process            0                   1                 2               3            4           5
rate               -1                  -1                -1              -1          -1          -1
------------------------------------------------------------------------------------
lumi               lnN    1.023        1.023             1.023           1.023      1.023       1.023
purewt             lnN    1.046        1.046             1.046           1.046      1.046       1.046
pdfrewt            lnN    -            1.010             1.010           1.010      1.010       1.010
topxsec            lnN    -            1.061             -               -          -           -
wxsec              lnN    -            -                 1.038           -          -           -
bc                 lnN    1.10         1.10              1.10            1.10       1.10        1.10
l                  lnN    1.07         1.07              1.07            1.07       1.07        1.07
mufr               lnN    1.12         1.12              1.12            1.12       1.12        1.12            
lep                lnN    1.03         1.03              1.03            1.03       1.03        1.03            
JES                lnN    1.024        1.024             1.024           1.024      1.024       1.024
JER                lnN    1.01         1.01              1.01            1.01       1.01        1.01            
'''

def MkDatacard(carddir, proc, mass, reg, ws):
    if not os.path.exists(os.path.join(os.getcwd(),carddir)):
        os.mkdir(carddir)

    if reg.startswith('1'):
        reg_name = reg.replace('1l', 'onel')
    elif reg.startswith('2'):
        reg_name = reg.replace('2l', 'twol')
    elif reg.startswith('3'):
        reg_name = reg.replace('3l', 'threel')
    else:
        reg_name = reg

    card = open(os.path.join(carddir, 'datacard_{0}_{1}_{2}.txt'.format(proc, mass, reg)), 'w')
    card_contents = re.sub('PROC', proc,    datacard_template)
    card_contents = re.sub('MASS', mass,    card_contents)
    card_contents = re.sub('REG_NAME',  reg_name,     card_contents)
    card_contents = re.sub('REG',  reg,     card_contents)
    card_contents = re.sub('WS',   ws,      card_contents)
  
    card.write(card_contents)
    card.close()


def GetProcs():
  
    chs = [
        'SMST2bW'  , 
        'SMST2tt'  ,  
        ]
    stop_mass = ['500', '600', '700']
    d=OrderedDict()
    for stop in stop_mass:
        d[stop] = []
    sigs=OrderedDict()
    for ch in chs:
        sigs[ch] = copy.deepcopy(d)
  
    sigs['SMST2bW']['500']        = [420, 440, 460, 480, 490]
    sigs['SMST2bW']['600']        = [590]
    sigs['SMST2bW']['700']        = [690]
  
    sigs['SMST2tt']['500']        = [420, 440, 460, 480, 490]
    sigs['SMST2tt']['600']        = [590]
    sigs['SMST2tt']['700']        = [690]
  
    print(sigs)
  
    lep_reg = ['1l']
    sub_reg = ['ge1ssv']
    #risr_reg = ['1p0', '0p98', '0p95', '0p9', '0p85', '0p8', 'le0p8']
    risr_reg = ['all']
     
    procs=OrderedDict([
        ('data_obs'            , OrderedDict()), 
        ('WJets2017'           , OrderedDict()), 
        ('TTJets2017'          , OrderedDict()), 
        ('ZJets2017'           , OrderedDict()), 
        ('VV2017'              , OrderedDict()), 
        ('DY2017'              , OrderedDict()), 
        ]) 
  
    for sig in sigs:
        for stop in sigs[sig]:
            lsp_masses = sigs[sig][stop]
            if len(lsp_masses) == 0: continue
            else:
                for lsp in lsp_masses:
                    procname = '{0}_{1}_{2}'.format(sig, stop, lsp)
                    procs[procname] = OrderedDict()
  
    for proc in procs:
        for lep in lep_reg:
            for sub in sub_reg:
                for risr in risr_reg:
                    if 'all' in risr:
                        procs[proc]['_'.join([lep, risr, sub])] = []
                    else:
                        procs[proc]['_'.join([lep, 'risr', risr, sub])] = []
  
    return procs, sigs


def MkWorkSpaceAndCards(fname, carddir):
  
    procs, sigs= GetProcs()
  
    print(procs)
  
    systs=[
        ]
  
    w = ROOT.RooWorkspace('w')
  
    #lo=0
    #hi=100
    #nbins=4
    lo=-0.5
    hi=34.5
    nbins=4
    bins = [0, 5, 10, 20, 40, 100]
    v_mass = ROOT.RooRealVar('m', 'm_{#perp} [GeV]', lo, hi)
    l_mass = ROOT.RooArgList(v_mass)
    s_mass = ROOT.RooArgSet(v_mass)
  
    fin = ROOT.TFile.Open(fname, 'READ')
  
    for proc, regs in procs.items():
        for reg in regs:
      
            print('>>>>>>>> MkWorkSpaceAndCards: Getting nominal histo for "{0}_MPE_{1}"'.format(proc, reg))
      
            h = fin.Get("{0}_MPE_{1}".format(proc, reg))
            #h.Rebin(2)
            if h.Integral() == 0:
                h.SetBinContent(1, 0.000001)

            nbins = h.GetNbinsX()
            if debug: print(">>>>>>>> MkWorkSpaceAndCards: h {0} has {1} bins".format(h.GetName(), nbins))

            lo = h.GetBinLowEdge(1)
            hi = h.GetBinLowEdge(nbins+1)
            v_mass.setBins(nbins)
        
            print('>>>>>>>> MkWorkSpaceAndCards: hname  {0}'.format(h.GetName()))
      
            roohist = ROOT.RooDataHist("{0}_MPE_{1}".format(proc, reg), "{0}_MPE_{1}".format(proc, reg), l_mass, h)
      
            if debug: print(">>>>>>>> MkWorkSpaceAndCards: roohist {0} has {1} bins".format(roohist.GetName(), roohist.numEntries()))
      
            getattr(w, 'import')(roohist)
      
            if 'Data' in proc: continue
      
            for syst in systs:
                print('>>>>>>>> MkWorkSpaceAndCards: Getting syst histo for proc {0} reg {1} syst {2}'.format(proc, reg, syst))
      
                hsystup = fin.Get("{0}_MPE_{1}_{2}Up".format(proc, reg, sys))
                hsystdown = fin.Get("{0}_MPE_{1}_{2}Down".format(proc, reg, sys))

                if not hsystup or not hsystdown: continue
      
                roohistup = ROOT.RooDataHist("{0}_MPE_{1}_{2}Up".format(proc, reg, syst), "{0}_MPE_{1}_{2}Up".format(proc, reg, sys), l_mass, hsystup)
                roohistdown = ROOT.RooDataHist("{0}_MPE_{1}_{2}Down".format(proc, reg, sys), "{0}_MPE_{1}_{2}Down".format(proc, reg, sys), l_mass, hsystdown)
      
                getattr(w, 'import')(roohistup)
                getattr(w, 'import')(roohistdown)
      
            ###Making datacard
            if re.match('SMS', proc) != None:
                print(proc.split('_'))
                mass = '_'.join(proc.split('_')[1:])
      
                print('Making datacard proc {0} mass {1} region {2}'.format(proc.split('_')[0], mass, reg))
                MkDatacard(carddir, proc.split('_')[0], mass, reg, fname.replace('template', 'workspace'))
    
    if not os.path.exists(carddir+'/workspaces'): os.mkdir(carddir+'/workspaces')
    w.writeToFile(os.path.join(carddir, fname.replace('template', 'workspace')))
  
    print("Workspace made and written to file!")


def CombineCards(carddir):
    print("Combining cards")
  
    procs, sigs = GetProcs()

    lep_reg = ['1l']
    sub_reg = ['ge1ssv']
    #risr_reg = ['1p0', '0p98', '0p95', '0p9', '0p85', '0p8', 'le0p8']
    risr_reg = ['all']
    regs = []

    for lep in lep_reg:
        for sub in sub_reg:
            for risr in risr_reg:
                if 'all' in risr:
                    regs.append('_'.join([lep, risr, sub]))
                else:
                    regs.append('_'.join([lep, 'risr', risr, sub]))
 
    pwd = os.getcwd()
    for proc in sigs: 
        stop_masses = sigs[proc]
        for stop in stop_masses:
            lsp_masses = sigs[proc][stop]
            for lsp in lsp_masses: 
                cards = ' '
                for reg in regs:
                    cards = cards + '{0}_{1}_{2}=datacard_{0}_{1}_{2}.txt '.format(proc, str(stop)+'_'+str(lsp), reg)
  
                cmd_template = 'combineCards.py  CARDS > datacard_PROC_STOP_LSP.txt'
                cmd_template = re.sub('CARDS', cards, cmd_template)
                cmd_template = re.sub('PROC', proc, cmd_template)
                cmd_template = re.sub('STOP', str(stop), cmd_template)
                cmd_template = re.sub('LSP', str(lsp), cmd_template)
                cmd_template = re.sub('DATACARDPATH', carddir, cmd_template)
                print(cmd_template)
  
                os.chdir(os.path.join(pwd,carddir))
                os.system(cmd_template)
                os.chdir(pwd) 
                os.system('text2workspace.py --channel-masks {0}/datacard_{1}_{2}.txt'.format(carddir, proc, stop+'_'+str(lsp)))
  
                submitAsymptotic(proc, stop, lsp, carddir)


import re, subprocess

asym_template='''#!/bin/bash


workdir=WORKDIR
rundir=RUNDIR

cd ${workdir}
eval `scramv1 runtime -sh`
cd ${rundir}

#combine -M FitDiagnostics -s 123456 -m STOP -d ../../CARDDIR/datacard_PROC_STOP_LSP.root --expectSignal 0 --cminFallbackAlgo Minuit,0.001 --rMin -1.0 --rMax 1.0 --robustFit=1 
combine -M AsymptoticLimits -s 123456 -m STOP -d ../../CARDDIR/datacard_PROC_STOP_LSP.root --run expected -v 5 --rMin -0.1 --rMax 5 --noFitAsimov >& AsymptoticLimits_PROC_STOP_LSP.log
mv higgsCombineTest.AsymptoticLimits.mHSTOP.root higgsCombineTest.AsymptoticLimits.PROC_STOP_LSP.root

#PostFitShapesFromWorkspace -m MASS\
#  -w ../../CARDDIR/datacard_PROC_STOP_LSP.root \
#  -o postfitshapes_PROC_STOP_LSP.root \
#  -f fitDiagnostics.root:fit_b \
#  --postfit \
#  --sampling \
#  --print \
#  >& postfit.log
cd ${workdir}
'''


def submitAsymptotic(sig, stop, lsp, carddir):

  #workdir = os.getcwd()
  workdir = os.getcwd()+'/lim_'+carddir
  if not os.path.exists(workdir): os.mkdir(workdir)
  rundir=os.path.join(workdir, 'limSig_{0}_{1}_{2}'.format(sig, stop, lsp))
  if not os.path.exists(rundir): os.mkdir(rundir)

  run_script = re.sub('WORKDIR', workdir, asym_template)
  run_script = re.sub('RUNDIR', rundir, run_script)
  run_script = re.sub('PROC', sig, run_script)
  run_script = re.sub('STOP', str(stop), run_script)
  run_script = re.sub('LSP', str(lsp), run_script)
  run_script = re.sub('CARDDIR', carddir, run_script)
  card = open(os.path.join(rundir, 'bash_{0}_{1}_{2}.sh'.format(sig, stop, lsp)), 'w')
  card.write(run_script)
  card.close()

  os.chdir(rundir)
  print('>>>>Running {}'.format(os.path.join(rundir, 'bash_{0}_{1}_{2}.sh'.format(sig, stop, lsp))))
  subprocess.check_output(['/bin/bash', os.path.join(rundir, 'bash_{0}_{1}_{2}.sh'.format(sig, stop, lsp))])

  os.chdir(workdir)


def main():

  usage='''
  Takes 2 arguments: (i) the selection templates, (ii) the dir to write the datacards.
  python MkWorkSpace.py templates/template.root datacards_selection/
  '''

  if len(sys.argv) < 2: 
    sys.exit(usage)

  import functools
  map(functools.partial(MkWorkSpaceAndCards, carddir=sys.argv[-1]), sys.argv[1:-1])
  CombineCards(sys.argv[-1])
  #runAsymptotic(sys.argv[-1])

if __name__ == "__main__":
  start = time.time() 
  main()
  stop = time.time()
  print('Run time: ', stop-start)
