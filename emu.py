import ROOT, math
import os
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine('.L Loader.C+') 
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetFillStyle(ROOT.kWhite)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetFrameBorderMode(ROOT.kWhite)
ROOT.gStyle.SetFrameFillColor(ROOT.kWhite)
ROOT.gStyle.SetCanvasBorderMode(ROOT.kWhite)
ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
ROOT.gStyle.SetPadBorderMode(ROOT.kWhite)
ROOT.gStyle.SetPadColor(ROOT.kWhite)
ROOT.gStyle.SetStatColor(ROOT.kWhite)
ROOT.gStyle.SetErrorX(0)

#First create the histos, then plot
#make_file is True if you want to create the histos, and then plot
#make_file is False if you want just to plot

make_file = True 


class MCSample:
    def __init__(self, name, crossSection, nEvents, color):
        self.name = name
        self.crossSection = crossSection
        self.color = color
        self.nEvents = nEvents
        self.effective_luminosity = self.nEvents/self.crossSection

class electron_object:
    def __init__(self, p4, charge):
        self.p4 = p4
        self.charge =charge

class muon_object:
    def __init__(self, p4, charge):
        self.p4 = p4
        self.charge = charge

class e_mu_object:
    def __init__(self, el, mu):
        self.el = el
        self.mu = mu
        self.p4 = el.p4 + mu.p4

conf = [ 'ss' , 'os' ] #same sign, opposite sign

hBase_mEE = {}
hBase_mEE['ss' ] = ROOT.TH1F('hBase_mEE_ss' , '', 200, 0, 2000)
hBase_mEE['os' ] = ROOT.TH1F('hBase_mEE_os', '',  200, 0, 2000)
for r in conf:
    hBase_mEE[r].GetXaxis().SetTitle('m(e#mu) [GeV]')
    hBase_mEE[r].Sumw2()

hBase_mEE['ss' ].GetYaxis().SetTitle('entries per 10 GeV' )
hBase_mEE['os'].GetYaxis().SetTitle('entries per 10 GeV')

samples = {}

samples['DYJetsToLL'    ] = MCSample('DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2'    ,  5.476e-5,  141842, ROOT.kMagenta  )
#samples['PHYS14_TT_20bx25'                  ] = MCSample('PHYS14_TT_20bx25'                  ,     689.1, 2423211, ROOT.kOrange)

if make_file:
    for sname in samples:
        s = samples[sname]
        s.hMEE = {}
        for r in conf:
            s.hMEE[r] = hBase_mEE[r].Clone('hMEE_%s_%s'%(sname,r))
            s.hMEE[r].SetFillColor(s.color)
            s.hMEE[r].SetLineWidth(0)

        if(sname=='DYJetsToLL'):
            s.tree= ROOT.TChain("IIHEAnalysis")
            for file in os.listdir("/user/aidan/public/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/"):
                if file.endswith(".root"):
                    s.tree.Add(str('/user/aidan/public/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/'+file))
        else:
            s.file = ROOT.TFile('the_other_path_for_tt')
            s.tree = s.file.Get('IIHEAnalysis')

        print s.name , s.tree.GetEntries()

        counter_entries=0

        for i in range(0,500000):#s.tree.GetEntries()):
            electrons= []
            muons    = []
            s.tree.GetEntry(i)

            if i%10000==0:
                print i , '/' , s.tree.GetEntries()

            if(sname=='DYJetsToLL'):# select only events with 2 taus, with mother Z
                counter_taus=0
                for i in range(0,s.tree.mc_n):
                    if(abs(s.tree.mc_pdgId[i])==15 and s.tree.mc_mother_pdgId[i][0]==23):
                        counter_taus=counter_taus + 1
                if (counter_taus <2):
                    continue

            ##ELECTRONS
            for j in range(0,s.tree.gsf_n):
                #Asking for HEEP selection;
                if s.tree.HEEP_cutflow51_total[j]==False:
                    continue
                el_Et  = s.tree.gsf_energy[j]/math.cosh(s.tree.gsf_eta[j])#*abs(math.sin(s.tree.gsf_theta[j])))
                el_eta = s.tree.gsf_eta[j]
                el_phi = s.tree.gsf_phi[j]
                el_E   = s.tree.gsf_energy[j]
                el_p4 = ROOT.TLorentzVector()
                el_p4.SetPtEtaPhiE(el_Et, el_eta, el_phi, el_E)
                el_charge = s.tree.gsf_charge[j]
                electrons.append(electron_object(el_p4,el_charge))
            # Sort by pt
            sorted(electrons, key=lambda e: e.p4.Pt())
            if len(electrons) == 0:
                continue

            ##MUONS
            for j in range(0,s.tree.mu_gt_n):#global tag muon
                ## Muon Selection
                selection=False
                selection= (s.tree.mu_gt_dxy[j]<2) and (s.tree.mu_gt_dz[j]<5)
                selection*=((s.tree.mu_gt_ptError[j]/s.tree.mu_gt_pt[j]) < 0.3)
                selection*= (s.tree.mu_numberOfMatchedStations[j] > 1) 

                #selection*=(s.tree.mu_numberOfValidTrackerHits[j] > 5) ##ALWAYS 0
                #selection*=(s.tree.mu_numberOfValidPixelHits[j] > 0) ##ALWAYS 0
                #(s.tree.mu_numberOfValidMuonHits[j]>0) ##ALWAYS 0
 
                if(not selection):
                    continue
                mu_Et  = abs(s.tree.mu_tevOptimized_pt[j]) 
                mu_eta = s.tree.mu_tevOptimized_eta[j]
                mu_phi = s.tree.mu_tevOptimized_phi[j]
                mu_E   = abs(s.tree.mu_tevOptimized_pt[j])*math.cosh(mu_eta) 
                mu_p4 = ROOT.TLorentzVector()
                mu_p4.SetPtEtaPhiE(mu_Et, mu_eta, mu_phi, mu_E)
                ##VETO IF MUON CLOSE TO ELE
                veto=False
                for ele in electrons:
                    dr=mu_p4.DeltaR(ele.p4)
                    if (dr<0.1 and ele.p4.Pt()>5):
                        veto=True
                if(veto):
                    continue

                mu_charge = s.tree.mu_tevOptimized_charge[j]
                muons.append(muon_object(mu_p4,mu_charge))
            # Sort by pt
            sorted(muons, key=lambda e: e.p4.Pt())
            if len(muons) ==0:
                continue
            counter_entries+=1
            E_mu = e_mu_object(electrons[0], muons[0])
            if(electrons[0].charge*muons[0].charge > 0):
                s.hMEE['ss'].Fill(E_mu.p4.M())
            else:
                s.hMEE['os'].Fill(E_mu.p4.M())

    print "number of entries with at least 1 muon and 1 ele", counter_entries
    file = ROOT.TFile('histograms_emu.root','RECREATE')

    ## Write the histos
    for sname in samples:
        s = samples[sname]
        for r in conf:
            s.hMEE[r].Write()
    file.Write()
    file.Close()
else: 
    file = ROOT.TFile('histograms_emu.root','READ')
    for sname in samples:
        s = samples[sname]
        s.hMEE = {}
        for r in conf:
            s.hMEE[r] = file.Get('hMEE_%s_%s'%(sname,r))
            s.hMEE[r].Scale(1000)
            s.hMEE[r].Scale(lumi/s.effective_luminosity)


############At this point you have rescaled your histos and you are ready to plot them#####################
canvas = ROOT.TCanvas('canvas','',100,100,800,600)
canvas.SetGridx()
canvas.SetGridy()

CMS_label_texts = {}
CMS_label_texts['normal'        ] = 'CMS'
CMS_label_texts['internal'      ] = 'CMS internal'
CMS_label_texts['workInProgress'] = 'CMS work in progress'
CMS_labels = {}
for t in CMS_label_texts:
    CMS_labels[t] = ROOT.TLatex(0.65, 0.945, CMS_label_texts[t])
    CMS_labels[t].SetNDC()
CMS_label = CMS_labels['internal']

lumi_label_texts = {}
lumi_label_texts['1'  ] = '#int L dt = 1 fb^{-1}'
lumi_label_texts['10'] = '#int L dt = 10 fb^{-1}'

lumi_labels = {}
for t in lumi_label_texts:
    lumi_labels[t] = ROOT.TLatex(0.5, 0.58, lumi_label_texts[t])
    lumi_labels[t].SetNDC()
lumi_label = lumi_labels['10']
lumi = 10 #for the 25 ns scenario, let's consider 10 fb-1

beam_label = ROOT.TLatex(0.25, 0.58, '#sqrt{s}=13 TeV')
beam_label.SetNDC()

x1 = 0.25
y1 = 0.85
x2 = 0.85
y2 = y1-0.2
legend = ROOT.TLegend(x1,y1,x2,y2)
legend.SetFillColor(0)
legend.SetBorderSize(0)
legend.SetShadowColor(0)
legend.SetNColumns(2)
legend.AddEntry(samples['DYJetsToLL'              ].hMEE[conf[0]], 'Z/#gamma#rightarrow#tau#tau', 'f')
#legend.AddEntry(samples['PHYS14_TT_20bx25'                  ].hMEE[conf[0]], "t#bar{t}"        , 'f')

hStack = {}

backgrounds = ['DYJetsToLL']

for r in conf:
    for l in ['lin','log']:
        hStack[r] = ROOT.THStack()
        for bname in backgrounds:
            hStack[r].Add(samples[bname].hMEE[r])
    
        min = 1e-3 if l=='log' else 0
        scale = 1e4 if l=='log' else 1.5
        hStack[r].SetMinimum(min)
        hStack[r].SetMaximum(scale*hStack[r].GetMaximum())
        hStack[r].Draw('hist')
        hStack[r].GetXaxis().SetTitle('m(e#mu) [GeV]')
        hStack[r].GetYaxis().SetTitle(samples[bname].hMEE[r].GetYaxis().GetTitle())
        hStack[r].GetYaxis().SetTitleOffset(1.25)
        hStack[r].Draw('hist')
        #for sname in samples: #Data
        #    if 'data' in sname:
        #        s = samples[sname]
        #        s.hMEE[r].SetFillColor(0)
        #        s.hMEE[r].SetLineWidth(2)
        #        s.hMEE[r].Draw('sames:hist')
        legend.Draw()
        CMS_label.Draw()
        lumi_label.Draw()
        beam_label.Draw()
        canvas.SetLogy(l=='log')
        canvas.Print('plots/M_emu_%s_%s.eps'%(r,l))


