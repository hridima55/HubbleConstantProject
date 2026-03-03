#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

#initial reading of CSV file with the distance data
import numpy as np

distances = np.genfromtxt('Distance_Mpc.csv', delimiter=',', dtype=None, encoding=None)

#discarding any data where the instrument response is not valid
value_to_remove = 'n'
columns_to_remove = np.any(distances == value_to_remove, axis=0) 
filtered_distances_array = distances[:, ~columns_to_remove] #this line of code was adapted from ChatGPT
print(filtered_distances_array)

velocities_array = [] #empty array, we will append each calculated value of redshifted velocity to the array to match to the distances array
velocity_uncertainties_array = []

#calculating the value of the redshifted value, v

#matching the observation number to the name of the spectral data .txt file
import os

directory_path = '/Users/hridimab/Downloads/Data(6)'
spectral_data_files = [f for f in os.listdir(directory_path) if f.endswith('.txt')] 
row_to_search = filtered_distances_array[0]
for i in row_to_search:
    file_name = f"{i}.txt"
    if file_name in spectral_data_files:
        file_path = os.path.join(directory_path, file_name)
        with open(file_path, 'r') as file: #reading the data
            spectral_data = file.readlines()[2:] #skipping the first two rows
            spectral_data_array = np.loadtxt(spectral_data) 
            
            frequency_values = spectral_data_array[:, 1]
            intensity_values = spectral_data_array[:, 2]
            
#initial guesses
            
            from scipy.optimize import curve_fit
            
            #defining the Gaussian function
            def gaussian(frequency_values, A, mu, sigma): #A = amplitude, mu = mean, sigma = standard deviation
                return A * np.exp(-(frequency_values - mu)**2 / (2 * sigma**2))
            
            #determining the initial guesses for the values of the slope and intercept of the linear part of the graph
            slope_guess, intercept_guess = np.polyfit(frequency_values, intensity_values, 1)
            
            #determining the initial guess for the value of the amplitude
            #removing the linear background to refine the estimate
            background = (slope_guess * frequency_values) + intercept_guess 
            intensity_no_background = intensity_values - background
            A_guess = np.median(intensity_no_background) #using median to estimate amplitude as the data is noisy
            
            #determining the intiial guess for the value of the mean
            mu_guess = frequency_values[np.argmax(intensity_no_background)]  #corresponding x-value of max y-value (i.e. the mean)
            from scipy.ndimage import gaussian_filter1d
            smoothed_intensity = gaussian_filter1d(intensity_no_background, sigma=2)
            #mu_guess = frequency_values[np.argmax(smoothed_intensity)]

            #determining the inital guess for the value of the standard deviation
            half_max = A_guess / 2
            #finding the full width at half maximum, FWHM, to guess a value for the standard deviation, sigma
            #here, the 'spread' of the Gaussian curve is being estimated by measuring the width at approximately half the maximum height (FWHM)
            #using the approximation: sigma = FWHM/2.355 (this idea and the following 9 lines of code were adapted from ChatGPT)
            half_max = np.max(smoothed_intensity) / 2
            x_half_max = frequency_values[smoothed_intensity > half_max]
            fwhm = x_half_max[-1] - x_half_max[0]
            sigma_guess = fwhm / 2.355
            #x_half_max = frequency_values[np.where(intensity_values > half_max)]
            #fwhm = x_half_max[-1] - x_half_max[0]  # Approximate width at half height
            #sigma_guess = fwhm / 2.355
            popt, pcov = curve_fit(gaussian, frequency_values, intensity_values, p0=[A_guess, mu_guess, sigma_guess],maxfev=100000)  #initial guesses for A, mu, sigma
            fitted_A, fitted_mu, fitted_sigma = popt #extract the fitted parameters
            print("Observation " + i + ":")
            print(" Initial Guesses: ")
            print(f"    Estimated Amplitude (A): {A_guess}")
            print(f"    Estimated Peak Location (mean): {mu_guess}")
            print(f"    Estimated Width (standard deviation): {sigma_guess}")
            print(f"    Estimated Slope: {slope_guess}")
            print(f"    Estimated Intercept: {intercept_guess}")
            
#fitting the data to Gaussian and linear function
           
            #the following 8 lines of code have been adapted from the solutions to the final advanced exercise in Core Worksheet 2
            def fit_function(frequency_values, A, mu, sigma, m, c):
               gaussian_part = A * np.exp(-(frequency_values - mu)**2 / (2 * sigma**2))
               linear_part = m * frequency_values + c
               return gaussian_part + linear_part
           
            initial_guesses = [A_guess, mu_guess, sigma_guess, slope_guess, intercept_guess] 
            
            po, pcov = curve_fit(fit_function, frequency_values, intensity_values, initial_guesses, maxfev=100000)

#plotting the initial guess line to see if it aligns with the data
            
            import matplotlib.pyplot as plt
            plt.plot(frequency_values, intensity_values, label="Data")
            plt.plot(
                frequency_values,
                fit_function(frequency_values, *initial_guesses),
                label="Initial Guess Model",
                color="orange",
                )  
            plt.xlabel("Frequency (Hz)")
            plt.ylabel("Spectra (A.U.)")
            plt.title("Initial Guess Model for Observation " + i)
            plt.legend()
            plt.show()

#plotting the graphs of frequency against intensity

            plt.plot(frequency_values, intensity_values, label="Spectral Data")
            plt.plot(frequency_values, fit_function(frequency_values, po[0], po[1], po[2], po[3], po[4]), label='Fit Results', color="orange")
            plt.xlabel("Frequency (Hz)")
            plt.ylabel("Spectra (A.U.)")
            plt.title("Graph for Observation " + i)
            plt.legend()
            
            plt.show()
            
#printing the value of the Gaussian mean (as this is equal to the observed frequency) and calculating the error
            gaussian_peak_index = np.argmax(intensity_values)
            gaussian_mean = frequency_values[gaussian_peak_index]
            uncertainty_in_mean = np.sqrt(pcov[1, 1])
            print (f" Observed Frequency: {gaussian_mean} +/- {uncertainty_in_mean} Hz")
            
#calculating the value of the observed wavelength using c = f_o * lambda_o
            observed_wavelength = (2.9979 * 10**8)/gaussian_mean
            uncertainty_in_wavelength = ((2.9979 * 10**8)/gaussian_mean**2) * uncertainty_in_mean
            print(f" Observed Wavelength: {observed_wavelength} +/- {uncertainty_in_wavelength} m")
    
#calculating the value of the redshifted velocity
            v = (2.9979 * 10**8)*((observed_wavelength)**2 - (486.1 * 10**-9)**2)/((observed_wavelength)**2 + (486.1 * 10**-9)**2)
            v = v/1000 #converting from m/s to km/s
            velocities_array.append(v)
            #calculating the uncertainty in velocity using the combining uncertainties formula
            uncertainty_in_v = (((2.9979 * 10**8)*((4*observed_wavelength*(486.1 * 10**-9)**2))))*(uncertainty_in_wavelength)/(observed_wavelength**2 + (486.1 * 10**-9)**2)
            uncertainty_in_v = uncertainty_in_v/1000 #converting from m/s to km/s
            velocity_uncertainties_array.append(uncertainty_in_v)
            print(f" Calculated Value of the Redshifted Velocity: {v} +/- {uncertainty_in_v} km/s")
            print(" ") #a line break to make the printed output look more tidy and easier to read
           
#matching each value of the redshifted velocity, v, to each value of distance, D
distance_values_array = filtered_distances_array[1] #isolates the distance values
distance_values_array = distance_values_array[1:]#removes the first header element

#plotting the final graph of v against D
import matplotlib.pyplot as plt

x = np.array(distance_values_array, dtype=float)
y = np.array(velocities_array, dtype=float)
yerr = np.array(velocity_uncertainties_array, dtype=float)
w = 1/(yerr**2)

x_weighted_mean = np.average(x, weights=w)
y_weighted_mean = np.average(y, weights=w)
slope = np.sum(w * (x - x_weighted_mean) * (y - y_weighted_mean)) / np.sum(w * (x - x_weighted_mean)**2)
intercept = y_weighted_mean - slope * x_weighted_mean

plt.errorbar(x, y, yerr=yerr, fmt='o', capsize=3)
plt.xlabel("Distance (Mpc)")
plt.ylabel("Redshifted Velocity (km/s)")
plt.title("Final Graph of Redshifted Velocity Against Distance")
line_of_best_fit = (slope * x) + intercept #plotting the line of best fit:
plt.plot(x, line_of_best_fit, color='orange')

plt.show()

#calculating the uncertainty in the value for H_0 using the covariance matrix
coefficients, covariance_matrix = np.polyfit(x, y, 1, w=w, cov=True)
uncertainty_in_H_0 = np.sqrt(covariance_matrix[0, 0])

#rounding the values of H_0 and the uncertainty to three significant figures
#the following function has been adapted from ChatGPT
import math

def round_to_significant_figures(number, sig_figs):
    if number == 0:
        return 0  #special case for 0
    else:
        return round(number, -int(math.floor(math.log10(abs(number)))) + (sig_figs - 1))

slope = round_to_significant_figures (slope, 3)
uncertainty_in_H_0 = round_to_significant_figures(uncertainty_in_H_0, 3)

print(f"The final value of Hubble's constant, H_0, is: {slope} +/- {uncertainty_in_H_0} (km/s)/Mpc")

           
            
                   


    
            


    







