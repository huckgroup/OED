from Model import ModelObject
"""A Super Main file meant to test a lot of individial models versus the pipeline that
was developed, a benchmark to spot inconsistancies essentially"""
class optimization_information:
    def __init__(self,theta = [],score = [],fitfunc = [],include = [],conditions = [],plot = []):
        self.theta      = theta
        self.score      = score
        self.fitfunc    = fitfunc
        self.include    = include
        self.conditions = conditions
        self.plot       = plot


def main():
    models = {}
    control = {}
    ##### you build a function that defines a space for a givenset of bounfaries
    description = "Phosphorylation model"
    
    """define parameters in model"""
    kinetic_parameters = ['k_cat13_AK_ADP',  'K_M13_ADP_AK','k_cat14_AK_AMP_ATP', 'K_M14_AMP_AK', 'K_M14_ATP_AK', 'k_cat15_PK_ADP_PEP', 'K_M15_ADP_PK', 'K_M15_PEP_PK', 'k_cat16_PK_PEP_UDP', 'K_M16_PEP_PK', 'K_M16_UDP_PK', 'k_cat17_PK_GDP_PEP', 'K_M17_GDP_PK', 'K_M17_PEP_PK', 'k_cat18_GMPK_ATP_GMP', 'K_M18_ATP_GMPK', 'K_M18_GMP_GMPK', 'k_cat19_UMPK_ATP_UMP', 'K_M19_ATP_UMPK', 'K_M19_UMP_UMPK', 'k_cat21_UPRT_PRPP_Ura'
                          ,'K_M21_PRPP_UPRT', 'K_M21_Ura_UPRT', 'k_cat22_APRT_Ade_PRPP', 'K_M22_Ade_APRT', 'K_M22_PRPP_APRT', 'k_cat26_PK_ATP_Pyr', 'K_M26_ATP_PK', 'K_M26_Pyr_PK', 'k_cat27_PK_Pyr_UTP', 'K_M27_Pyr_PK', 'K_M27_UTP_PK', 'k_cat28_PK_GTP_Pyr', 'K_M28_GTP_PK', 'K_M28_Pyr_PK', 'k_cat29_GMPK_ADP_GDP', 'K_M29_ADP_GMPK', 'K_M29_GDP_GMPK', 'k_cat30_UMPK_ADP_UDP', 'K_M30_ADP_UMPK', 'K_M30_UDP_UMPK']
   
    enzymes = ['PK', 'UMPK', 'UPRT', 'GMPK', 'APRT', 'AK']
    for i in enzymes:
        kinetic_parameters.append(i)
    
    """define control parameters in model"""
    control_parameters = ['Ade_in', 'Gua_in', 'PEP_in', 'PRPP_in', 'Ura_in', 'ATP_in']

    """give an initial guess for the value of each parameter, set to 1 here"""
    fixed = {k: 1 for k in kinetic_parameters} 
    for i in enzymes:
        fixed[i] = 1
    fixed.update({c: 10 for c in control_parameters})
    for i in ['ADP_in', 'AMP_in', 'ATP_in', 'Ade_in', 'GDP_in', 'GMP_in', 'GTP_in', 'Gua_in', 'PEP_in', 'PPi_in', 'PRPP_in', 'Pyr_in', 'UDP_in', 'UMP_in', 'UTP_in', 'Ura_in']:
        fixed[i] = 0
    fixed.update({"kf": 0.125})

    """define the upper and lower boundary of each parameter"""
    boundaries = {'Ade_in':(0.01,10),'Gua_in':(0.01,1),'PEP_in':(0.01,10),'PRPP_in':(0.1,10),'Ura_in':(0.01,1),'ATP_in':(0.01,1)}
    boundaries.update({k: (0.001, 50) for k,v in fixed.items() if k not in control_parameters})
    boundaries.update({"kf": (0.01, 0.25)})

    """define equations of you system"""
    stringmodel = ''' -2*ADP**2*AK*k_cat13_AK_ADP/(K_M13_ADP_AK*(ADP/K_M13_ADP_AK |+| 1)) - ADP*GDP*GMPK*k_cat29_GMPK_ADP_GDP/(K_M29_ADP_GMPK*K_M29_GDP_GMPK*(ADP/K_M29_ADP_GMPK |+| 1)*(GDP/K_M29_GDP_GMPK |+| 1)) - ADP*kf - ADP*UDP*UMPK*k_cat30_UMPK_ADP_UDP/(K_M30_ADP_UMPK*K_M30_UDP_UMPK*(1 |+| UDP/K_M30_UDP_UMPK)*(ADP/K_M30_ADP_UMPK |+| 1)) - ADP*PEP*PK*k_cat15_PK_ADP_PEP/(K_M15_ADP_PK*K_M15_PEP_PK*(1 |+| PEP/K_M15_PEP_PK)*(ADP/K_M15_ADP_PK |+| 1)) + ADP_in*kf + 2*AK*AMP*ATP*k_cat14_AK_AMP_ATP/(K_M14_AMP_AK*K_M14_ATP_AK*(AMP/K_M14_AMP_AK |+| 1)*(ATP/K_M14_ATP_AK |+| 1)) + ATP*GMP*GMPK*k_cat18_GMPK_ATP_GMP/(K_M18_ATP_GMPK*K_M18_GMP_GMPK*(ATP/K_M18_ATP_GMPK |+| 1)*(GMP/K_M18_GMP_GMPK |+| 1)) + ATP*PK*Pyr*k_cat26_PK_ATP_Pyr/(K_M26_ATP_PK*K_M26_Pyr_PK*(1 |+| Pyr/K_M26_Pyr_PK)*(ATP/K_M26_ATP_PK |+| 1)) + ATP*UMP*UMPK*k_cat19_UMPK_ATP_UMP/(K_M19_ATP_UMPK*K_M19_UMP_UMPK*(1 |+| UMP/K_M19_UMP_UMPK)*(ATP/K_M19_ATP_UMPK |+| 1)) 
    +ADP**2*AK*k_cat13_AK_ADP/(K_M13_ADP_AK*(ADP/K_M13_ADP_AK |+| 1)) - AK*AMP*ATP*k_cat14_AK_AMP_ATP/(K_M14_AMP_AK*K_M14_ATP_AK*(AMP/K_M14_AMP_AK |+| 1)*(ATP/K_M14_ATP_AK |+| 1)) - AMP*kf + AMP_in*kf + APRT*Ade*PRPP*k_cat22_APRT_Ade_PRPP/(K_M22_Ade_APRT*K_M22_PRPP_APRT*(1 |+| PRPP/K_M22_PRPP_APRT)*(Ade/K_M22_Ade_APRT |+| 1)) 
    +ADP**2*AK*k_cat13_AK_ADP/(K_M13_ADP_AK*(ADP/K_M13_ADP_AK |+| 1)) + ADP*GDP*GMPK*k_cat29_GMPK_ADP_GDP/(K_M29_ADP_GMPK*K_M29_GDP_GMPK*(ADP/K_M29_ADP_GMPK |+| 1)*(GDP/K_M29_GDP_GMPK |+| 1)) + ADP*UDP*UMPK*k_cat30_UMPK_ADP_UDP/(K_M30_ADP_UMPK*K_M30_UDP_UMPK*(1 |+| UDP/K_M30_UDP_UMPK)*(ADP/K_M30_ADP_UMPK |+| 1)) + ADP*PEP*PK*k_cat15_PK_ADP_PEP/(K_M15_ADP_PK*K_M15_PEP_PK*(1 |+| PEP/K_M15_PEP_PK)*(ADP/K_M15_ADP_PK |+| 1)) - AK*AMP*ATP*k_cat14_AK_AMP_ATP/(K_M14_AMP_AK*K_M14_ATP_AK*(AMP/K_M14_AMP_AK |+| 1)*(ATP/K_M14_ATP_AK |+| 1)) - ATP*GMP*GMPK*k_cat18_GMPK_ATP_GMP/(K_M18_ATP_GMPK*K_M18_GMP_GMPK*(ATP/K_M18_ATP_GMPK |+| 1)*(GMP/K_M18_GMP_GMPK |+| 1)) - ATP*kf - ATP*PK*Pyr*k_cat26_PK_ATP_Pyr/(K_M26_ATP_PK*K_M26_Pyr_PK*(1 |+| Pyr/K_M26_Pyr_PK)*(ATP/K_M26_ATP_PK |+| 1)) - ATP*UMP*UMPK*k_cat19_UMPK_ATP_UMP/(K_M19_ATP_UMPK*K_M19_UMP_UMPK*(1 |+| UMP/K_M19_UMP_UMPK)*(ATP/K_M19_ATP_UMPK |+| 1)) + ATP_in*kf 
    -APRT*Ade*PRPP*k_cat22_APRT_Ade_PRPP/(K_M22_Ade_APRT*K_M22_PRPP_APRT*(1 |+| PRPP/K_M22_PRPP_APRT)*(Ade/K_M22_Ade_APRT |+| 1)) - Ade*kf + Ade_in*kf 
    -ADP*GDP*GMPK*k_cat29_GMPK_ADP_GDP/(K_M29_ADP_GMPK*K_M29_GDP_GMPK*(ADP/K_M29_ADP_GMPK |+| 1)*(GDP/K_M29_GDP_GMPK |+| 1)) + ATP*GMP*GMPK*k_cat18_GMPK_ATP_GMP/(K_M18_ATP_GMPK*K_M18_GMP_GMPK*(ATP/K_M18_ATP_GMPK |+| 1)*(GMP/K_M18_GMP_GMPK |+| 1)) - GDP*kf - GDP*PEP*PK*k_cat17_PK_GDP_PEP/(K_M17_GDP_PK*K_M17_PEP_PK*(1 |+| PEP/K_M17_PEP_PK)*(GDP/K_M17_GDP_PK |+| 1)) + GDP_in*kf + GTP*PK*Pyr*k_cat28_PK_GTP_Pyr/(K_M28_GTP_PK*K_M28_Pyr_PK*(1 |+| Pyr/K_M28_Pyr_PK)*(GTP/K_M28_GTP_PK |+| 1)) 
    +ADP*GDP*GMPK*k_cat29_GMPK_ADP_GDP/(K_M29_ADP_GMPK*K_M29_GDP_GMPK*(ADP/K_M29_ADP_GMPK |+| 1)*(GDP/K_M29_GDP_GMPK |+| 1)) - ATP*GMP*GMPK*k_cat18_GMPK_ATP_GMP/(K_M18_ATP_GMPK*K_M18_GMP_GMPK*(ATP/K_M18_ATP_GMPK |+| 1)*(GMP/K_M18_GMP_GMPK |+| 1)) - GMP*kf + GMP_in*kf 
    +GDP*PEP*PK*k_cat17_PK_GDP_PEP/(K_M17_GDP_PK*K_M17_PEP_PK*(1 |+| PEP/K_M17_PEP_PK)*(GDP/K_M17_GDP_PK |+| 1)) - GTP*kf - GTP*kf - GTP*PK*Pyr*k_cat28_PK_GTP_Pyr/(K_M28_GTP_PK*K_M28_Pyr_PK*(1 |+| Pyr/K_M28_Pyr_PK)*(GTP/K_M28_GTP_PK |+| 1)) + GTP_in*kf 
    -ADP*PEP*PK*k_cat15_PK_ADP_PEP/(K_M15_ADP_PK*K_M15_PEP_PK*(1 |+| PEP/K_M15_PEP_PK)*(ADP/K_M15_ADP_PK |+| 1)) + ATP*PK*Pyr*k_cat26_PK_ATP_Pyr/(K_M26_ATP_PK*K_M26_Pyr_PK*(1 |+| Pyr/K_M26_Pyr_PK)*(ATP/K_M26_ATP_PK |+| 1)) - GDP*PEP*PK*k_cat17_PK_GDP_PEP/(K_M17_GDP_PK*K_M17_PEP_PK*(1 |+| PEP/K_M17_PEP_PK)*(GDP/K_M17_GDP_PK |+| 1)) + GTP*PK*Pyr*k_cat28_PK_GTP_Pyr/(K_M28_GTP_PK*K_M28_Pyr_PK*(1 |+| Pyr/K_M28_Pyr_PK)*(GTP/K_M28_GTP_PK |+| 1)) - PEP*kf + PEP_in*kf + PK*Pyr*UTP*k_cat27_PK_Pyr_UTP/(K_M27_Pyr_PK*K_M27_UTP_PK*(1 |+| Pyr/K_M27_Pyr_PK)*(1 |+| UTP/K_M27_UTP_PK)) - PEP*PK*UDP*k_cat16_PK_PEP_UDP/(K_M16_PEP_PK*K_M16_UDP_PK*(1 |+| PEP/K_M16_PEP_PK)*(1 |+| UDP/K_M16_UDP_PK)) 
    +APRT*Ade*PRPP*k_cat22_APRT_Ade_PRPP/(K_M22_Ade_APRT*K_M22_PRPP_APRT*(1 |+| PRPP/K_M22_PRPP_APRT)*(Ade/K_M22_Ade_APRT |+| 1)) - PPi*kf + PPi_in*kf + PRPP*UPRT*Ura*k_cat21_UPRT_PRPP_Ura/(K_M21_PRPP_UPRT*K_M21_Ura_UPRT*(1 |+| PRPP/K_M21_PRPP_UPRT)*(1 |+| Ura/K_M21_Ura_UPRT)) 
    -APRT*Ade*PRPP*k_cat22_APRT_Ade_PRPP/(K_M22_Ade_APRT*K_M22_PRPP_APRT*(1 |+| PRPP/K_M22_PRPP_APRT)*(Ade/K_M22_Ade_APRT |+| 1)) - PRPP*kf + PRPP_in*kf - PRPP*UPRT*Ura*k_cat21_UPRT_PRPP_Ura/(K_M21_PRPP_UPRT*K_M21_Ura_UPRT*(1 |+| PRPP/K_M21_PRPP_UPRT)*(1 |+| Ura/K_M21_Ura_UPRT)) 
    +ADP*PEP*PK*k_cat15_PK_ADP_PEP/(K_M15_ADP_PK*K_M15_PEP_PK*(1 |+| PEP/K_M15_PEP_PK)*(ADP/K_M15_ADP_PK |+| 1)) - ATP*PK*Pyr*k_cat26_PK_ATP_Pyr/(K_M26_ATP_PK*K_M26_Pyr_PK*(1 |+| Pyr/K_M26_Pyr_PK)*(ATP/K_M26_ATP_PK |+| 1)) + GDP*PEP*PK*k_cat17_PK_GDP_PEP/(K_M17_GDP_PK*K_M17_PEP_PK*(1 |+| PEP/K_M17_PEP_PK)*(GDP/K_M17_GDP_PK |+| 1)) - GTP*PK*Pyr*k_cat28_PK_GTP_Pyr/(K_M28_GTP_PK*K_M28_Pyr_PK*(1 |+| Pyr/K_M28_Pyr_PK)*(GTP/K_M28_GTP_PK |+| 1)) - Pyr*kf + Pyr_in*kf - PK*Pyr*UTP*k_cat27_PK_Pyr_UTP/(K_M27_Pyr_PK*K_M27_UTP_PK*(1 |+| Pyr/K_M27_Pyr_PK)*(1 |+| UTP/K_M27_UTP_PK)) + PEP*PK*UDP*k_cat16_PK_PEP_UDP/(K_M16_PEP_PK*K_M16_UDP_PK*(1 |+| PEP/K_M16_PEP_PK)*(1 |+| UDP/K_M16_UDP_PK)) 
    -ADP*UDP*UMPK*k_cat30_UMPK_ADP_UDP/(K_M30_ADP_UMPK*K_M30_UDP_UMPK*(1 |+| UDP/K_M30_UDP_UMPK)*(ADP/K_M30_ADP_UMPK |+| 1)) + ATP*UMP*UMPK*k_cat19_UMPK_ATP_UMP/(K_M19_ATP_UMPK*K_M19_UMP_UMPK*(1 |+| UMP/K_M19_UMP_UMPK)*(ATP/K_M19_ATP_UMPK |+| 1)) - UDP*kf + UDP_in*kf + PK*Pyr*UTP*k_cat27_PK_Pyr_UTP/(K_M27_Pyr_PK*K_M27_UTP_PK*(1 |+| Pyr/K_M27_Pyr_PK)*(1 |+| UTP/K_M27_UTP_PK)) - PEP*PK*UDP*k_cat16_PK_PEP_UDP/(K_M16_PEP_PK*K_M16_UDP_PK*(1 |+| PEP/K_M16_PEP_PK)*(1 |+| UDP/K_M16_UDP_PK))
    +ADP*UDP*UMPK*k_cat30_UMPK_ADP_UDP/(K_M30_ADP_UMPK*K_M30_UDP_UMPK*(1 |+| UDP/K_M30_UDP_UMPK)*(ADP/K_M30_ADP_UMPK |+| 1)) - ATP*UMP*UMPK*k_cat19_UMPK_ATP_UMP/(K_M19_ATP_UMPK*K_M19_UMP_UMPK*(1 |+| UMP/K_M19_UMP_UMPK)*(ATP/K_M19_ATP_UMPK |+| 1)) - UMP*kf + UMP_in*kf + PRPP*UPRT*Ura*k_cat21_UPRT_PRPP_Ura/(K_M21_PRPP_UPRT*K_M21_Ura_UPRT*(1 |+| PRPP/K_M21_PRPP_UPRT)*(1 |+| Ura/K_M21_Ura_UPRT)) 
    -UTP*kf + UTP_in*kf - PK*Pyr*UTP*k_cat27_PK_Pyr_UTP/(K_M27_Pyr_PK*K_M27_UTP_PK*(1 |+| Pyr/K_M27_Pyr_PK)*(1 |+| UTP/K_M27_UTP_PK)) + PEP*PK*UDP*k_cat16_PK_PEP_UDP/(K_M16_PEP_PK*K_M16_UDP_PK*(1 |+| PEP/K_M16_PEP_PK)*(1 |+| UDP/K_M16_UDP_PK)) 
    -Ura*kf + Ura_in*kf - PRPP*UPRT*Ura*k_cat21_UPRT_PRPP_Ura/(K_M21_PRPP_UPRT*K_M21_Ura_UPRT*(1 |+| PRPP/K_M21_PRPP_UPRT)*(1 |+| Ura/K_M21_Ura_UPRT)) '''

    symbols = ['-','+',')','(','*','**','|+|',"|-|","/",'\n']
    reconstructed = [('| + |','|+|'),('| - |','|-|'),("*  *","**")]
    for i in symbols:
        stringmodel = stringmodel.replace(i," " + i + " " )
    for i,j in reconstructed:
        stringmodel = stringmodel.replace(i,j)
    
    print(stringmodel)
    
    #define states of equations, make sure they are in the same order (or that you define them seperately in a dictionary with equations)
    states      = ['ADP', 'AMP', 'ATP', 'Ade', 'GDP', 'GMP', 'GTP', 'PEP', 'PPi', 'PRPP', 'Pyr', 'UDP', 'UMP', 'UTP', 'Ura']
    observables = ['ADP', 'ATP', 'Ade', 'GDP', 'GMP', 'GTP', 'UDP',  'UTP']
    
    #define the scoring function you wish to use
    conditions = {e: 0.1 for e in enzymes}
    conditions.update({c: 10 for c in control_parameters})
    conditions.update({"kf": 0.125})

    #the modelname is simultenously the name of the folder where all the information abou the model will be stored
    modelname             = 'NSP_model'
    include               = []
    models[len(models)]   = ModelObject(stringmodel,states,boundaries,fixed,observables = observables,name = modelname,control_parameters = control_parameters)
    control[len(control)] = optimization_information(conditions = conditions)
    
    return models,control
