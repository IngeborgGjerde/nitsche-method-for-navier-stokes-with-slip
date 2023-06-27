from fenics import *
import numpy as np
from pathlib import Path
def write_convergence_table(error, file_name, variable_name, parameters, gammas, meshes, hs, hgammas):

    desc = '_beta'+str(parameters['beta'])

    caption = 'with $\\beta=%1.1f$ and '%parameters['beta'] # table caption
    if parameters['cyl_refinement']==True:
        desc += '_refined'
        caption += 'refined mesh'
    else:
        desc += '_uniform'
        caption += 'uniform mesh'

    if parameters['beta']==-2.0:
        desc += '_potflow'
        caption += ', potential flow'
    else:
        desc += '_mansol'
        caption += ', manufactured solution'

    if parameters['normal_choice']=='projected_normal':
        desc += '_projnormal'
        caption += ', projected normal'
    else:
        desc += '_discnormal'
        caption += ', discrete normal'

    file_path = Path(file_name+desc)
    file_path.parent.mkdir(parents=True, exist_ok=True)    
    ffile = open(file_name+desc, "w+")

    # Make preamble for table
    table_preamble = '\\begin{table} \n \\begin{center} \n' + \
                          '\\caption{'  + variable_name + caption+ '} \n'

    # Write table preamble and column headings etc
    ffile.write(table_preamble)
    ffile.write('\\begin{tabular}{|l|*{%g}{l}|}\hline'%(len(gammas)) + '\n')
    ffile.write('\\backslashbox{$h_\Omega$, \, $h_\Gamma$}{$\gamma$}' + '\n ')

    # Make columns for each value of gamma
    gamma_line = ''
    for gamma in gammas:
        temp = '&\\makebox[3em]{%.1e}'%gamma
        gamma_line += temp.replace('e-0', 'e-').replace('e+00', '')
    ffile.write(gamma_line + '\\\ \hline\hline \n')

    # Write down errors at each iteration, recording conv. rate as we go 
    convs = np.zeros((len(meshes)-1, len(gammas)))

    if parameters['cyl_refinement']: h=hgammas
    else: h=hs

    for ix_N in range(len(meshes)):
        strr = ''
        # Write column for mesh size
        strr += "%.3f"%hs[ix_N] + ", %.3f"%hgammas[ix_N]

        # Iterate through columns for gamma, writing down error and calculating the convergence rate
        for ix_gamma in range(len(gammas)):

            temp = '&  ' + "%.1e"%error[ix_N][ix_gamma]
            strr += temp.replace('e-0', 'e-').replace('e+00', '')

            if ix_N>0:
                error_dx = np.log(error[ix_N][ix_gamma])-np.log(error[ix_N-1][ix_gamma])
                h_dx = np.log(h[ix_N])-np.log(h[ix_N-1])
                convs[ix_N-1][ix_gamma] = error_dx/h_dx

                strr += (' (' +  '%.1f'%convs[ix_N-1][ix_gamma] + ')')

        strr += ' \\\ \n'
        ffile.write(strr)

    # Write convergence rates
    conv_avgs = np.average(convs, axis=0)
    strr = '\hline \n  avg. rate '
    for ix_gamma, conv_rate in enumerate(conv_avgs):
        strr += '&  %.1f '%conv_avgs[ix_gamma]
    strr += '\n \\\ \hline \n' 
    ffile.write(strr)

    # Write down end of table
    ffile.write('\end{tabular} \n \end{center} \n \end{table}')
    ffile.close()