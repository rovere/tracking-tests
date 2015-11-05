from ROOT import TFile, gDirectory, TH1F, TH2F
import math

# open the file
inputFile = TFile( 'ntuple.root' )

# retrieve the ntuple of interest
trkTree = gDirectory.Get( 'trkTree/tree' )
entries = trkTree.GetEntriesFast()

for jentry in xrange( entries ):
    
    # get the next tree in the chain and verify
    ientry = trkTree.LoadTree( jentry )
    if ientry < 0: break
    # copy next entry into memory and verify
    nb = trkTree.GetEntry( jentry )
    if nb <= 0: continue

    # loop over simulated tracks
    for itp in range(0,len(trkTree.sim_pt)):
        q =  trkTree.sim_q[itp]
        px =  trkTree.sim_px[itp]
        py =  trkTree.sim_py[itp]
        pz =  trkTree.sim_pz[itp]
        pt = trkTree.sim_pt[itp]
        eta = trkTree.sim_eta[itp]
        phi = trkTree.sim_phi[itp]
        prodx =  trkTree.sim_prodx[itp]
        prody =  trkTree.sim_prody[itp]
        prodz =  trkTree.sim_prodz[itp]
        if math.fabs(eta)>0.5: continue
        if math.fabs(prodz)>0.1: continue        
        trkIdx = trkTree.sim_trkIdx[itp]
        trkAlgo = -1
        if trkIdx >= 0: trkAlgo = trkTree.trk_algo[trkIdx]
        hits = []
        #sim hits
        for ipix in trkTree.sim_pixelIdx[itp]:
            x = trkTree.pix_xsim[ipix]
            y = trkTree.pix_ysim[ipix]
            z = trkTree.pix_zsim[ipix]
            r = math.sqrt(x*x + y*y)
            eta = -math.log( math.tan( math.atan2(r,z)/2) )
            phi = math.atan2(y,x)
            det = 2-trkTree.pix_isBarrel[ipix]
            lay = trkTree.pix_lay[ipix]
            hits.append( (det,lay,x,y,z,r,eta,phi) )
            #print ( "\tpixel det=%i lay=%i sim pos=%f,%f,%f r=%f eta=%f phi=%f" % (det,lay,x,y,z,r,eta,phi) )
        for istr in trkTree.sim_stripIdx[itp]:
            #if trkTree.str_isStereo[istr]==1: continue
            x = trkTree.str_xsim[istr]
            y = trkTree.str_ysim[istr]
            z = trkTree.str_zsim[istr]
            r = math.sqrt(x*x + y*y)
            eta = -math.log( math.tan( math.atan2(r,z)/2) )
            phi = math.atan2(y,x)
            det = trkTree.str_det[istr]
            lay = trkTree.str_lay[istr]
            hits.append( (det,lay,x,y,z,r,eta,phi) )
            #print ( "\tstrip det=%i lay=%i sim pos=%f,%f,%f r=%f eta=%f phi=%f" % (det,lay,x,y,z,r,eta,phi) )
        prev_det = -1
        prev_lay = -1
        totLayers = 0
        for hit in sorted(hits, key=lambda hit: hit[5]):
            if hit[0]==prev_det and hit[1]==prev_lay: continue
            prev_det = hit[0]
            prev_lay = hit[1]
            totLayers = totLayers+1
        if totLayers < 13: continue
        #now we can print the info
        print ( "simTrack %f %f %f %f %f %f %i pt=%f trkIdx=%i algo=%i" % (prodx,prody,prodz,px,py,pz,q,pt,trkIdx,trkAlgo) )
        prev_det = -1
        prev_lay = -1
        for hit in sorted(hits, key=lambda hit: hit[5]):
            if hit[0]==prev_det and hit[1]==prev_lay: continue
            print  ( "simHit %f %f %f %f %f" % (hit[2],hit[3],hit[4],hit[5],hit[6]) )
            prev_det = hit[0]
            prev_lay = hit[1]
