#Sample usage:

import rdkit, rdkit.Chem
from rdkit.Chem.PandasTools import LoadSDF
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
import openpyxl
PandasTools.RenderImagesInAllDataFrames(images=True)
from rdkit.Chem import Recap,BRICS
import rdkit.Chem.AllChem as AllChem

import sklearn
import xgboost as xgb


smiles_file = pd.read_csv('input_file.csv')

voting_clf_fp = pickle.load(open("bms_ligands_maximum_voting_clf.sav", 'rb'))
nBits = 2048
ecfp6_name = [f'Bit_{i}' for i in range(nBits)]
remover = SaltRemover()

smiles = [x.strip().replace('*','') for x in list(smiles_file['SMILES'])]

prediction_dict = {'smiles':[], 'mol':[], 'prediction_class_ml':[],'predicted_activity_ml':[],\
                   'predicted_pIC50_ridge_regression':[], 'qed':[]}

for index, row in predicted_actives_best.iterrows():
    
    try:
        smi = row['SMILES']
        
        mol = Chem.MolFromSmiles(smi)
        mol = remover(Chem.AddHs(mol))

        fplist = []
        fp = AllChem.GetMorganFingerprintAsBitVect(mol,4,nBits=nBits)
        fplist.append(fp)
        ecfp6_bits = [list(l) for l in fplist]
        ecfp6_df = pd.DataFrame(ecfp6_bits, columns=ecfp6_name)

        prediction_class = voting_clf_fp.predict(ecfp6_df)[0]

        if prediction_class == 0:
            predicted_activity = 'Inactive'
        elif prediction_class == 1:
            predicted_activity = 'Active'

        del ecfp6_df, fplist, fp, ecfp6_bits

        prediction_dict['smiles'].append(smi)
        prediction_dict['mol'].append(mol)
        prediction_dict['prediction_class_ml'].append(prediction_class)
        prediction_dict['predicted_activity_ml'].append(predicted_activity)
        prediction_dict['qed'].append(row['qed'])
        
    except:
        pass

prediction_df = pd.DataFrame.from_dict(prediction_dict, orient = 'columns')

prediction_df[prediction_df['prediction_class_ml'] == 1].drop(columns = ['mol']).to_excel('predictions_output_actives.xlsx')
prediction_df.head()
