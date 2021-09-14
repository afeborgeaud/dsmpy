import numpy as np
from Goes_attenuation_models import anelastic_properties, Q7g, fn_pressure
import matplotlib.pyplot as plt
import os
"""
Choose the filename, and define the name of the output file
"""
data_directory = "./data_directory"
result_directory = "./result_directory"

file_number = [37]
areal_number = [1,128,256,384,512,640,758]

for a in range(0,len(file_number)):
    for b in range(0,len(areal_number)):
        prefix_name = "perplex_binary_2d1e8_R104"
        file_number_string = "%05.0f"%file_number[a]
        areal_number_string = "%05.0f"%areal_number[b]
        file_name = f"{prefix_name}_{areal_number_string}_{file_number_string}"
        
        input_file_name = os.path.join(data_directory,file_name+".dat")
        output_file_name = os.path.join(result_directory,file_name+".out")

        print("input_file_name:",input_file_name)
        print("output_file_name:",output_file_name)

        data = np.loadtxt(input_file_name)     

        # Here we define the parameters for the attenuation model
        frequency = 1.
        Q_model = Q7g
        dT_Q_constant_above_solidus=0

        # Set up required datat arrays
        pressure = np.zeros(len(data))
        anelastic_Vp = np.zeros(len(data))
        anelastic_Vs = np.zeros(len(data))
        Q_S = np.zeros(len(data))
        Q_K = np.zeros(len(data))
        Q_P = np.zeros(len(data))
        T_solidus = np.zeros(len(data))

        # Process the elastic input using the anlelasticity model


        for i,d in enumerate(data):
            depth, temperature, c, rho, elastic_Vs, vb, elastic_Vp = d
            pressure[i] = fn_pressure(depth*1000.)

            dict = anelastic_properties(elastic_Vp, elastic_Vs, depth*1000., temperature,
                                        frequency, Q_model, dT_Q_constant_above_solidus)
            anelastic_Vp[i] = dict.V_P
            anelastic_Vs[i] = dict.V_S
            Q_S[i] = dict.Q_S
            Q_K[i] = dict.Q_K
            Q_P[i] = dict.Q_P
            T_solidus[i] = dict.T_solidus

        # Create a single 2D data array
        depth, temperature, c, rho, elastic_Vs, vb, elastic_Vp = data.T
        output = np.array([depth, temperature, c, rho*1000., elastic_Vp, elastic_Vs,
                        anelastic_Vp, anelastic_Vs, Q_S, Q_K, Q_P, T_solidus]).T

        # Save the file, including a header with units
        header = ('depth (km), T (K), c, rho (kg/m^3), elastic_Vp (km/s), elastic_Vs (km/s),'
        ' anelastic_Vp (km/s), anelastic_Vs (km/s), Q_S, Q_K, Q_P, T_solidus (K)')
        np.savetxt(output_file_name, output, header=header, fmt='%.5e')

        if True:
            fig = plt.figure()
            ax = [fig.add_subplot(2, 2, i) for i in range(1, 5)]

            ax[0].plot(depth, pressure/1.e9)
            ax[0].set_ylabel('Pressure (GPa)')
            ax[0].set_xlabel('Depth (km)')


            ax[1].plot(depth, temperature)
            ax[1].set_ylabel('Temperature (K)')
            ax[1].set_xlabel('Depth (km)')

            ax[2].plot(depth, Q_S, label='$Q_S$')
            ax[2].plot(depth, Q_K, label='$Q_K$')
            ax[2].plot(depth, Q_P, label='$Q_P$')
            ax[2].set_yscale('log')

            ax[2].set_ylabel('$Q$')
            ax[2].set_xlabel('Depth (km)')
            ax[2].legend()

            ax[3].plot(depth, c)
            ax[3].set_ylabel('Composition')
            ax[3].set_xlabel('Depth (km)')

            fig.set_tight_layout(True)
            fig.savefig(os.path.join(result_directory,file_name+".pdf"), dpi = 300, bbox_inches="tight")
            plt.show()
