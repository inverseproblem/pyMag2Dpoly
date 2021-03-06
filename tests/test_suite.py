
import sys
import numpy as np

sys.path.append("./")
import mag2dpoly as mag

########################################################

def test1():

    Jind = mag.MagnetizVector(mod=4.9,Ideg=90.0,Ddeg=45.0)

    Jrem = mag.MagnetizVector(mod=3.1,Ideg=45.0,Ddeg=0.0)

    northxax = 90.0

    Nobs = 101
    xzobs = np.transpose(np.vstack(( np.linspace(0.0,100.0,Nobs), -1.0*np.ones(Nobs))))

    vertices = np.array([ [35.0, 50.0],
                          [65.0, 50.0],
                          [80.0, 35.0],
                          [65.0, 20.0],
                          [35.0, 20.0],
                          [20.0, 35.0] ])


    ## array of arrays
    nbod = 1
    bodyindices = np.empty(shape=(nbod,), dtype=np.object)
    inds = range(6)
    bodyindices[0] = np.array(inds)

    pbody = mag.MagPolyBodies2D(bodyindices,vertices)

    ## compute total field (superposition of tmag of bodies)
    forwardtype = ["talwani","talwani_red","krav","wonbev"]
    nfwd = len(forwardtype)
    tmag = np.zeros((xzobs.shape[0],nfwd))
    nbody = pbody.bo.size

    Jindv = np.array([Jind]) 
    Jremv = np.array([Jrem])

    for (f,fwdt) in enumerate(forwardtype):
        tmag[:,f] = mag.tmagpolybodies2Dgen(xzobs,Jindv,Jremv,northxax,pbody,forwardtype[f])

    ##=========================================

    tmag_ref = np.array([ -116.16336912423309,
                          -107.60475622349567,
                          -98.17723975959588,
                          -87.82662842366796,
                          -76.49775256995866,	
                          -64.13484643609536,
                          -50.682001858897294,
                          -36.08369908036598,
                          -20.285419949994587,
                          -3.2343484367466506,
                          15.119837140202215,
                          34.82407638693401,
                          55.92094064483901,
                          78.44754640867708,
                          102.43436783270268,
                          127.90394975470824,
                          154.8695254168912,
                          183.3335469464661,
                          213.28614208527696,
                          244.70351807578447,
                          277.5463434505119,
                          311.758151040835,
                          347.2638208530303,
                          383.9682191188879,
                          421.7550886836475,
                          460.4863039069605,
                          500.0016173506611,
                          540.119031663933,
                          580.6359235938513,
                          621.331023429952,
                          661.9673092115332,
                          702.2958101588549,
                          742.0602315222982,
                          781.0022216070385,
                          818.8670137044077,
                          855.4091061413785,
                          890.3976071295564,
                          923.62087771655,
                          954.890158626987,
                          984.0419590104054,
                          1010.9391031766961,
                          1035.4704568268496,
                          1057.5494684492162,
                          1077.1117497867892,
                          1094.111973515768,
                          1108.5203855942902,
                          1120.3192191143874,
                          1129.499264418154,
                          1136.0568061804884,
                          1139.991090454424,
                          1141.302439274996,
                          1139.9910904544242,
                          1136.0568061804884,
                          1129.499264418154,
                          1120.3192191143867,
                          1108.52038559429,
                          1094.111973515767,
                          1077.1117497867892,
                          1057.5494684492162,
                          1035.4704568268494,
                          1010.9391031766954,
                          984.0419590104062,
                          954.8901586269874,
                          923.6208777165492,
                          890.3976071295564,
                          855.4091061413781,
                          818.8670137044072,
                          781.0022216070377,
                          742.0602315222977,
                          702.295810158854,
                          661.9673092115327,
                          621.3310234299516,
                          580.6359235938512,
                          540.1190316639322,
                          500.0016173506604,
                          460.48630390696064,
                          421.75508868364693,
                          383.9682191188882,
                          347.2638208530297,
                          311.7581510408341,
                          277.5463434505112,
                          244.70351807578382,
                          213.28614208527665,
                          183.33354694646584,
                          154.86952541689124,
                          127.90394975470798,
                          102.43436783270306,
                          78.4475464086772,
                          55.92094064483902,
                          34.82407638693425,
                          15.119837140202181,
                          -3.234348436746581,
                          -20.285419949994342,
                          -36.083699080365705,
                          -50.68200185889743,
                          -64.13484643609462,
                          -76.49775256995939,
                          -87.82662842366803,
                          -98.17723975959623,
                          -107.60475622349514,
                          -116.16336912423309])

    ##=========================================

    same = [np.allclose(tmag[:,i],tmag_ref) for i in range(tmag.shape[1])]

    print("all(same): ",all(same))

    if all(same) :
        return True
    else :
        for i in range(nfwd):
            if same[i]==False : 
                print("Failed for ",forwardtype[i])
            
        return False

################################################################
################################################################

def test2():

    Jind = mag.MagnetizVector(mod=4.9,Ideg=90.0,Ddeg=45.0)

    Jrem = mag.MagnetizVector(mod=3.1,Ideg=45.0,Ddeg=0.0)

    northxax = 90.0

    Nobs = 101
    xzobs = np.transpose(np.vstack(( np.linspace(0.0,100.0,Nobs), -1.0*np.ones(Nobs))))

    vertices = np.array([[39.9767752377201, 13.278991954030298],
                         [43.21168888652307, 13.402291449170173],
                         [43.317652312473136, 14.739290679698811],  
                         [43.546649812203974, 13.203638738574346],  
                         [43.91560484766585, 13.172384962659983],  
                         [44.14344907863333, 13.366977100309308],  
                         [44.40792260208563, 14.433477629153117],  
                         [44.69197601158024, 13.41736322165884],   
                         [43.051984435135324, 12.28972077523347],   
                         [41.64196146829342, 11.82578388310931],   
                         [40.0, 10.0],                
                         [51.593762919577486, 12.985398228772771], 
                         [52.893379913843816, 11.074145035188602],
                         [53.03238448156978, 10.749216912218833],
                         [53.17394928686575, 10.785849813929147],  
                         [53.21679397938986, 10.936364994289892],  
                         [53.31684510201022, 10.901007653562868],  
                         [52.164999382907524, 9.566098324208415],
                         [51.69760101144281, 10.0],                
                         [60.559336241924306, 25.07742274155419],
                         [65.61740252955528, 25.405022500973224], 
                         [66.6522482157798, 24.330905878795054],
                         [66.99370462700534, 22.430809598897174],
                         [67.95964200423342, 20.871100638551827], 
                         [68.41821282760033, 20.46729028637734],  
                         [67.67999306568265, 20.187645182621033], 
                         [67.42746927845577, 19.369846341974842], 
                         [65.71882888409377, 17.507867891105946], 
                         [65.63379439530597, 16.85561419292242],  
                         [65.52497931177504, 14.899927178489119], 
                         [64.20104498482468, 12.937540781866952],
                         [62.56031603648637, 11.598110231234838], 
                         [60.78820427692065, 10.0]])


    ## array of arrays
    nbod = 3
    bodyindices = np.empty(shape=(nbod,), dtype=np.object)
    ind1 = range(11)
    ind2 = range(11,19)
    ind3 = range(19,33)
    bodyindices[0] = np.array(ind1)
    bodyindices[1] = np.array(ind2)
    bodyindices[2] = np.array(ind3)

    pbody = mag.MagPolyBodies2D(bodyindices,vertices)

    ## compute total field (superposition of tmag of bodies)
 
    ## compute total field (superposition of tmag of bodies)
    forwardtype = ["talwani","talwani_red","krav","wonbev"]
    nfwd = len(forwardtype)
    tmag = np.zeros((xzobs.shape[0],nfwd))
    nbody = pbody.bo.size

    for (f,fwdt) in enumerate(forwardtype):
        for i in range(nbody):
            tmag[:,f] += mag.tmagpoly2Dgen(xzobs,Jind,Jrem,northxax,pbody.bo[i],forwardtype[f])

    ###########################################################################

    tmag_ref = [-25.18101203031008,
                -25.814111679386286,
                -26.46591188468871,
                -27.136513821069393,
                -27.82589071680914,
                -28.533854477346786,
                -29.260014742480834,
                -30.003728645158226,
                -30.764039135123586,
                -31.53959923300231,
                -32.32857897143186,
                -33.12855103994676,
                -33.93635025992559,
                -34.74790095972702,
                -35.558005092427436,
                -36.360082555206354,
                -37.14585368518767,
                -37.90495244621838,
                -38.624457626684254,
                -39.28832887474771,
                -39.87673534453047,
                -40.36526832574021,
                -40.724037390920856,
                -40.91666522359578,
                -40.89922361225488,
                -40.61919790412695,
                -40.01463677488,
                -39.013746213150505,
                -37.53532637130842,
                -35.49062356397606,
                -32.78735204013638,
                -29.336765353626436,
                -25.06459486138753,
                -19.926213598986948,
                -13.925267157081024,
                -7.133055532339174,
                0.2966799431282361,
                8.125321786345921,
                16.047468856320766,
                23.73409183413679,
                30.896898575646777,
                37.35947087475763,
                43.11364135233236,
                48.34131017064244,
                53.39281394546982,
                58.7274836881353,
                64.83295011470159,
                72.14348149663033,
                80.97541229042082,
                91.49257447124981,
                103.70790485858134,
                117.51814892303605,
                132.75663244957155,
                149.23869197241677,
                166.77366986321647,
                185.1314145984579,
                203.9747839694488,
                222.7894273624335,
                240.84727102364297,
                257.2295606565289,
                270.9152383168951,
                280.91830730494206,
                286.440275178761,
                286.99698347216184,
                282.4869194996864,
                273.1878553816839,
                259.6912319962327,
                242.79908349000044,
                223.41189102240213,
                202.4297310575936,
                180.67896175091312,
                158.86732337391078,
                137.5640994379945,
                117.19911113444167,
                98.07391651937287,
                80.37955931028789,
                64.21668987444241,
                49.615331309029074,
                36.552740278605995,
                24.96865481639862,
                14.7777648040935,
                5.879555320167113,
                -1.834170357150894,
                -8.473724523855406,
                -14.147563163141468,
                -18.959577242278417,
                -23.00736030015876,
                -26.381212717810055,
                -29.163688024210938,
                -31.429526675334152,
                -33.245857032715975,
                -34.67257173728606,
                -35.76281076807948,
                -36.5635008447628,
                -37.11591518972078,
                -37.45622870164328,
                -37.616051936809455,
                -37.6229334959912,
                -37.50082493720261,
                -37.27050556150594,
                -36.949966659975665]

    ##=================================================
  
    same = [np.allclose(tmag[:,i],tmag_ref) for i in range(tmag.shape[1])]

    print("all(same): ",all(same))

    if all(same) :
        return True
    else :
        for i in range(nfwd):
            if same[i]==False : 
                print("Failed for ",forwardtype[i])
        return False

#############################################################

if __name__=="__main__" :

    test1()
    test2()
