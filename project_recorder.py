import yaml
import sys
import read_function






class inDrop_project:
    yaml=''
    sample_index={}
    input_dir=''
    output_dir=''
    file_name={'Cellbarcode_1':'',
    'Cellbarcode_2':'',
    'library_index':'',
    'Biological_read':''}
    def __init__(self,project_yaml):
        self.yaml=yaml.load(project_yaml)
        self.sample_index=self.yaml['sample_index']
        self.input_dir=self.yaml['input_dir']
        self.output_dir=self.yaml['output_dir']
        Lanepresence=all(['[Lane]' not in name for name in [self.yaml['Cellbarcode_1'],self.yaml['Cellbarcode_2'],self.yaml['library_index'],self.yaml['Biological_read']]])
        if self.yaml['Lane'] is not None and Lanepresence is True:
            self.file_name['Cellbarcode_1']=self.yaml['Cellbarcode_1']
            self.file_name['Cellbarcode_2']=self.yaml['Cellbarcode_2']
            self.file_name['library_index']=self.yaml['library_index']
            self.file_name['Biological_read']=self.yaml['Biological_read']
        else:
            if self.yaml['Lane'] is None and Lanepresence is True:
                sys.exit('The file name input indicate that there are multiple lanes, but the Lane information is not provided.')
            else:
                for key in self.file_name:
                    self.file_name[key]=[]
                for lane in self.yaml['Lane']:
                    Cellbarcode_1=self.yaml['Cellbarcode_1'].split('[Lane]')
                    Cellbarcode_2=self.yaml['Cellbarcode_2'].split('[Lane]')
                    library_index=self.yaml['library_index'].split('[Lane]')
                    Biological_read=self.yaml['Biological_read'].split('[Lane]')
                    self.file_name['Cellbarcode_1'].append('%s%s%s'%(Cellbarcode_1[0],lane,Cellbarcode_1[1]))
                    self.file_name['Cellbarcode_2'].append('%s%s%s'%(Cellbarcode_2[0],lane,Cellbarcode_2[1]))
                    self.file_name['library_index'].append('%s%s%s'%(library_index[0],lane,library_index[1]))
                    self.file_name['Biological_read'].append('%s%s%s'%(Biological_read[0],lane,Biological_read[1]))


