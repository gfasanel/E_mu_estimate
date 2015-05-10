def fix_HEEP_cuts(tree, metatree,index):
    eta = abs(tree.gsf_eta[index])
    dEtaIn    = abs(tree.HEEP_cutflow41_dEtaIn_value[index])
    dEtaIn_51 = abs(tree.HEEP_cutflow51_dEtaIn_value[index])
    HOverE    = abs(tree.HEEP_cutflow41_HOverE_value[index])
    ET = tree.HEEP_cutflow41_Et_value[index]
    E = tree.gsf_energy[index]
    
    tree.HEEP_cutflow41_dEtaIn[index]      = True
    tree.HEEP_cutflow50_50ns_dEtaIn[index] = True
    tree.HEEP_cutflow50_25ns_dEtaIn[index] = True
    tree.HEEP_cutflow50_dEtaIn[index]      = True
    tree.HEEP_cutflow51_dEtaIn[index]      = True
    
    tree.HEEP_cutflow41_HOverE[index]      = True
    tree.HEEP_cutflow50_50ns_HOverE[index] = True
    tree.HEEP_cutflow50_25ns_HOverE[index] = True
    tree.HEEP_cutflow50_HOverE[index]      = True
    tree.HEEP_cutflow51_HOverE[index]      = True
    
    if eta < metatree.HEEP_cutflow_barrelEtaUpper_41:
        dEtaIn_threshold_41 = 0.005
        HOverE_threshold_41 = 0.05
        if dEtaIn > dEtaIn_threshold_41:
            tree.HEEP_cutflow41_dEtaIn[index] = False
        if HOverE > HOverE_threshold_41:
            tree.HEEP_cutflow41_HOverE[index] = False
    elif eta > metatree.HEEP_cutflow_endcapEtaLower_41 and eta < metatree.HEEP_cutflow_endcapEtaUpper_41:
        dEtaIn_threshold_41 = 0.005
        HOverE_threshold_41 = 0.05
        if dEtaIn > dEtaIn_threshold_41:
            tree.HEEP_cutflow41_dEtaIn[index] = False
        dEtaIn_threshold_41 = 0.005
        if dEtaIn > HOverE_threshold_41:
            tree.HEEP_cutflow41_HOverE[index] = False
    
    if eta < metatree.HEEP_cutflow_barrelEtaUpper_50:
        dEtaIn_threshold_50_50ns = max( 0.016-1e-4*ET , 0.004 )
        dEtaIn_threshold_50_25ns = max( 0.016-1e-4*ET , 0.004 )
        dEtaIn_threshold_50      = max( 0.016-1e-4*ET , 0.004 )
        dEtaIn_threshold_51 = 0.004
        
        HOverE_threshold_50_50ns = 2.0/E + 0.05
        HOverE_threshold_50_25ns = 2.0/E + 0.05
        HOverE_threshold_50      = 2.0/E + 0.05
        HOverE_threshold_51      = 2.0/E + 0.05
        
        if dEtaIn > dEtaIn_threshold_50_50ns:
            tree.HEEP_cutflow50_50ns_dEtaIn[index] = False
        if dEtaIn > dEtaIn_threshold_50_25ns:
            tree.HEEP_cutflow50_25ns_dEtaIn[index] = False
        if dEtaIn > dEtaIn_threshold_50:
            tree.HEEP_cutflow50_dEtaIn[index]      = False
        if dEtaIn_51 > dEtaIn_threshold_51:
            tree.HEEP_cutflow51_dEtaIn[index]      = False
        
        if HOverE > HOverE_threshold_50_50ns:
            tree.HEEP_cutflow50_50ns_HOverE[index] = False
        if HOverE > HOverE_threshold_50_25ns:
            tree.HEEP_cutflow50_25ns_HOverE[index] = False
        if HOverE > HOverE_threshold_50:
            tree.HEEP_cutflow50_HOverE[index]      = False
        if HOverE > HOverE_threshold_51:
            tree.HEEP_cutflow51_HOverE[index]      = False
            
    elif eta > metatree.HEEP_cutflow_endcapEtaLower_50 and eta < metatree.HEEP_cutflow_endcapEtaUpper_50:
        dEtaIn_threshold_50_50ns = 0.02
        dEtaIn_threshold_50_25ns = max( 0.015-8.5e-5*ET , 0.006 )
        dEtaIn_threshold_50      = max( 0.015-8.5e-5*ET , 0.006 )
        dEtaIn_threshold_51 = 0.006
        
        HOverE_threshold_50_50ns = 12.5/E + 0.05
        HOverE_threshold_50_25ns = 12.5/E + 0.05
        HOverE_threshold_50      = 12.5/E + 0.05
        HOverE_threshold_51      = 12.5/E + 0.05
        
        if dEtaIn > dEtaIn_threshold_50_50ns:
            tree.HEEP_cutflow50_50ns_dEtaIn[index] = False
        if dEtaIn > dEtaIn_threshold_50_25ns:
            tree.HEEP_cutflow50_25ns_dEtaIn[index] = False
        if dEtaIn > dEtaIn_threshold_50:
            tree.HEEP_cutflow50_dEtaIn[index]      = False
        if dEtaIn_51 > dEtaIn_threshold_50_25ns:
            tree.HEEP_cutflow51_dEtaIn[index]      = False
        
        if HOverE > HOverE_threshold_50_50ns:
            tree.HEEP_cutflow50_50ns_HOverE[index] = False
        if HOverE > HOverE_threshold_50_25ns:
            tree.HEEP_cutflow50_25ns_HOverE[index] = False
        if HOverE > HOverE_threshold_50:
            tree.HEEP_cutflow50_HOverE[index]      = False
        if HOverE > HOverE_threshold_51:
            tree.HEEP_cutflow51_HOverE[index]      = False
    
    tree.HEEP_cutflow41_ID[index] = True
    if tree.HEEP_cutflow41_EcalDriven[index]==False:
        tree.HEEP_cutflow41_ID[index] = False
    if tree.HEEP_cutflow41_dEtaIn[index]==False:
        tree.HEEP_cutflow41_ID[index] = False
    if tree.HEEP_cutflow41_dPhiIn[index]==False:
        tree.HEEP_cutflow41_ID[index] = False
    if tree.HEEP_cutflow41_HOverE[index]==False:
        tree.HEEP_cutflow41_ID[index] = False
    if tree.HEEP_cutflow41_SigmaIetaIeta[index]==False:
        tree.HEEP_cutflow41_ID[index] = False
    if tree.HEEP_cutflow41_E1x5OverE5x5[index]==False:
        tree.HEEP_cutflow41_ID[index] = False
    if tree.HEEP_cutflow41_E2x5OverE5x5[index]==False:
        tree.HEEP_cutflow41_ID[index] = False
    if tree.HEEP_cutflow41_missingHits[index]==False:
        tree.HEEP_cutflow41_ID[index] = False
    if tree.HEEP_cutflow41_dxyFirstPV[index]==False:
        tree.HEEP_cutflow41_ID[index] = False
    
    tree.HEEP_cutflow41_total[index] = tree.HEEP_cutflow41_acceptance[index] and tree.HEEP_cutflow41_ID[index] and tree.HEEP_cutflow41_isolation[index]
    
    
    
    tree.HEEP_cutflow50_50ns_ID[index] = True
    if tree.HEEP_cutflow50_50ns_EcalDriven[index]==False:
        tree.HEEP_cutflow50_50ns_ID[index] = False
    if tree.HEEP_cutflow50_50ns_dEtaIn[index]==False:
        tree.HEEP_cutflow50_50ns_ID[index] = False
    if tree.HEEP_cutflow50_50ns_dPhiIn[index]==False:
        tree.HEEP_cutflow50_50ns_ID[index] = False
    if tree.HEEP_cutflow50_50ns_HOverE[index]==False:
        tree.HEEP_cutflow50_50ns_ID[index] = False
    if tree.HEEP_cutflow50_50ns_SigmaIetaIeta[index]==False:
        tree.HEEP_cutflow50_50ns_ID[index] = False
    if tree.HEEP_cutflow50_50ns_E1x5OverE5x5[index]==False:
        tree.HEEP_cutflow50_50ns_ID[index] = False
    if tree.HEEP_cutflow50_50ns_E2x5OverE5x5[index]==False:
        tree.HEEP_cutflow50_50ns_ID[index] = False
    if tree.HEEP_cutflow50_50ns_missingHits[index]==False:
        tree.HEEP_cutflow50_50ns_ID[index] = False
    if tree.HEEP_cutflow50_50ns_dxyFirstPV[index]==False:
        tree.HEEP_cutflow50_50ns_ID[index] = False
    
    tree.HEEP_cutflow50_50ns_total[index] = tree.HEEP_cutflow50_50ns_acceptance[index] and tree.HEEP_cutflow50_50ns_ID[index] and tree.HEEP_cutflow50_50ns_isolation[index]
    
    
    
    tree.HEEP_cutflow50_25ns_ID[index] = True
    if tree.HEEP_cutflow50_25ns_EcalDriven[index]==False:
        tree.HEEP_cutflow50_25ns_ID[index] = False
    if tree.HEEP_cutflow50_25ns_dEtaIn[index]==False:
        tree.HEEP_cutflow50_25ns_ID[index] = False
    if tree.HEEP_cutflow50_25ns_dPhiIn[index]==False:
        tree.HEEP_cutflow50_25ns_ID[index] = False
    if tree.HEEP_cutflow50_25ns_HOverE[index]==False:
        tree.HEEP_cutflow50_25ns_ID[index] = False
    if tree.HEEP_cutflow50_25ns_SigmaIetaIeta[index]==False:
        tree.HEEP_cutflow50_25ns_ID[index] = False
    if tree.HEEP_cutflow50_25ns_E1x5OverE5x5[index]==False:
        tree.HEEP_cutflow50_25ns_ID[index] = False
    if tree.HEEP_cutflow50_25ns_E2x5OverE5x5[index]==False:
        tree.HEEP_cutflow50_25ns_ID[index] = False
    if tree.HEEP_cutflow50_25ns_missingHits[index]==False:
        tree.HEEP_cutflow50_25ns_ID[index] = False
    if tree.HEEP_cutflow50_25ns_dxyFirstPV[index]==False:
        tree.HEEP_cutflow50_25ns_ID[index] = False
    
    tree.HEEP_cutflow50_25ns_total[index] = tree.HEEP_cutflow50_25ns_acceptance[index] and tree.HEEP_cutflow50_25ns_ID[index] and tree.HEEP_cutflow50_25ns_isolation[index]
    
    
    
    tree.HEEP_cutflow50_ID[index] = True
    if tree.HEEP_cutflow50_EcalDriven[index]==False:
        tree.HEEP_cutflow50_ID[index] = False
    if tree.HEEP_cutflow50_dEtaIn[index]==False:
        tree.HEEP_cutflow50_ID[index] = False
    if tree.HEEP_cutflow50_dPhiIn[index]==False:
        tree.HEEP_cutflow50_ID[index] = False
    if tree.HEEP_cutflow50_HOverE[index]==False:
        tree.HEEP_cutflow50_ID[index] = False
    if tree.HEEP_cutflow50_SigmaIetaIeta[index]==False:
        tree.HEEP_cutflow50_ID[index] = False
    if tree.HEEP_cutflow50_E1x5OverE5x5[index]==False:
        tree.HEEP_cutflow50_ID[index] = False
    if tree.HEEP_cutflow50_E2x5OverE5x5[index]==False:
        tree.HEEP_cutflow50_ID[index] = False
    if tree.HEEP_cutflow50_missingHits[index]==False:
        tree.HEEP_cutflow50_ID[index] = False
    if tree.HEEP_cutflow50_dxyFirstPV[index]==False:
        tree.HEEP_cutflow50_ID[index] = False
        
    tree.HEEP_cutflow50_total[index] = tree.HEEP_cutflow50_acceptance[index] and tree.HEEP_cutflow50_ID[index] and tree.HEEP_cutflow50_isolation[index]
    
    
    
    tree.HEEP_cutflow51_ID[index] = True
    if tree.HEEP_cutflow51_EcalDriven[index]==False:
        tree.HEEP_cutflow51_ID[index] = False
    if tree.HEEP_cutflow51_dEtaIn[index]==False:
        tree.HEEP_cutflow51_ID[index] = False
    if tree.HEEP_cutflow51_dPhiIn[index]==False:
        tree.HEEP_cutflow51_ID[index] = False
    if tree.HEEP_cutflow51_HOverE[index]==False:
        tree.HEEP_cutflow51_ID[index] = False
    if tree.HEEP_cutflow51_SigmaIetaIeta[index]==False:
        tree.HEEP_cutflow51_ID[index] = False
    if tree.HEEP_cutflow51_E1x5OverE5x5[index]==False:
        tree.HEEP_cutflow51_ID[index] = False
    if tree.HEEP_cutflow51_E2x5OverE5x5[index]==False:
        tree.HEEP_cutflow51_ID[index] = False
    if tree.HEEP_cutflow51_missingHits[index]==False:
        tree.HEEP_cutflow51_ID[index] = False
    if tree.HEEP_cutflow51_dxyFirstPV[index]==False:
        tree.HEEP_cutflow51_ID[index] = False
        
    tree.HEEP_cutflow51_total[index] = tree.HEEP_cutflow51_acceptance[index] and tree.HEEP_cutflow51_ID[index] and tree.HEEP_cutflow51_isolation[index]
