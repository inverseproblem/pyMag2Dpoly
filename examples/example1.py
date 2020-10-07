
import sys
import numpy as np

sys.path.append("../")
import mag2dpoly as mag


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
tmag = np.zeros(xzobs.shape[0])
nbody = pbody.bo.size
    
forwardtype = "talwani"
    
for i in range(nbody):
    tmag += mag.tmagpoly2Dgen(xzobs,Jind,Jrem,northxax,pbody.bo[i],forwardtype)

    

##################################

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

##################################

assert np.allclose(tmag,tmag_ref)
print("np.allclose(tmag,tmag_ref): ",np.allclose(tmag,tmag_ref))
