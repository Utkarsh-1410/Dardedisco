import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import speed_of_light as c
from scipy.signal import windows
from mpl_toolkits.mplot3d import Axes3D

def simulate_full_ajisai():
    # S-band Radar parameters
    fc = 3e9  # 3 GHz
    wavelength = c / fc
    desired_resolution = 0.7  # 70 cm
    bw = c / (2 * desired_resolution)  # 214 MHz
    prf = 500  # Hz
    n_pulses = 256  # Increased for better angular resolution
    pulse_width = 20e-6  # 20 μs for longer range coverage
    sample_rate = 2 * bw  # 428 MHz
    n_samples = int(pulse_width * sample_rate)
    
    # AJISAI parameters
    satellite_diameter = 2.15  # m
    n_reflectors = 318  # Actual number of corner reflectors
    rotation_rate = 0.1  # rad/sec
    
    # Generate all 318 corner reflector positions (Golden Spiral distribution)
    indices = np.arange(0, n_reflectors, dtype=float) + 0.5
    theta = np.arccos(1 - 2*indices/n_reflectors)  # Polar angle
    phi = np.pi * (1 + 5**0.5) * indices  # Azimuthal angle
    
    # Convert to Cartesian coordinates
    x = satellite_diameter/2 * np.sin(theta) * np.cos(phi)
    y = satellite_diameter/2 * np.sin(theta) * np.sin(phi)
    z = satellite_diameter/2 * np.cos(theta)
    
    # Visualize the satellite geometry
    fig = plt.figure(figsize=(18, 6))
    
    # 3D View
    ax1 = fig.add_subplot(131, projection='3d')
    ax1.scatter(x, y, z, s=20, c='r')
    ax1.set_title('AJISAI Corner Reflector Distribution\n(318 Triangular Reflectors)')
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    ax1.set_zlabel('Z (m)')
    
    # Generate radar signal
    t = np.arange(n_samples) / sample_rate
    slope = bw / pulse_width
    tx_sig = np.exp(1j * np.pi * slope * t**2) * windows.tukey(n_samples, alpha=0.25)
    
    # Initialize echo matrix
    echo_matrix = np.zeros((n_samples, n_pulses), dtype=complex)
    
    # Radar viewing geometry (slant range 1500 km)
    radar_range = 1500e3  # 1500 km
    
    for pulse_idx in range(n_pulses):
        rot_angle = rotation_rate * pulse_idx / prf
        
        for i in range(n_reflectors):
            # Rotate the reflector
            x_rot = x[i] * np.cos(rot_angle) - y[i] * np.sin(rot_angle)
            y_rot = x[i] * np.sin(rot_angle) + y[i] * np.cos(rot_angle)
            
            # Calculate range to reflector (including z-axis variation)
            range_to_reflector = radar_range + z[i]  # Simplified model
            range_delay = 2 * range_to_reflector / c
            
            # Doppler from rotation
            doppler = 4 * np.pi * y_rot * rotation_rate / wavelength
            
            # Create echo
            phase = 2 * np.pi * fc * range_delay + doppler
            echo = np.sqrt(10) * tx_sig * np.exp(1j * phase)  # Assume 10 m² RCS
            
            # Add to echo matrix
            idx = int(np.round(range_delay * sample_rate))
            if idx < n_samples:
                end_idx = min(idx + len(echo), n_samples)
                echo_matrix[idx:end_idx, pulse_idx] += echo[:end_idx-idx]
    
    # Apply windowing
    range_window = windows.taylor(n_samples, nbar=5, sll=30)
    doppler_window = windows.hamming(n_pulses)
    echo_matrix = echo_matrix * range_window[:, np.newaxis] * doppler_window[np.newaxis, :]
    
    # ISAR processing
    range_fft = np.fft.fft(echo_matrix, axis=0)
    doppler_fft = np.fft.fftshift(np.fft.fft(range_fft, axis=1), axes=1)
    
    # Create ISAR image
    isar_image = 20 * np.log10(np.abs(doppler_fft) + 1e-10)
    isar_image = isar_image - np.max(isar_image)  # Normalize
    
    # ISAR Image Plot
    ax2 = fig.add_subplot(132)
    ax2.imshow(isar_image,
               extent=[-prf/2, prf/2, 0, n_samples/sample_rate*c/2],
               aspect='auto',
               cmap='jet',
               vmin=-30, vmax=0)
    ax2.set_title('S-band ISAR Image\n(3 GHz, 70 cm Resolution)')
    ax2.set_xlabel('Doppler Frequency (Hz)')
    ax2.set_ylabel('Range (m)')
    
    # Zoomed ISAR View
    ax3 = fig.add_subplot(133)
    ax3.imshow(isar_image,
               extent=[-prf/2, prf/2, radar_range-5, radar_range+5],  # Zoom to ±5m
               aspect='auto',
               cmap='jet',
               vmin=-30, vmax=0)
    ax3.set_title('Zoomed View: Satellite Structure')
    ax3.set_xlabel('Doppler Frequency (Hz)')
    ax3.set_ylabel('Range (m)')
    
    plt.tight_layout()
    plt.show()

# Run the simulation
simulate_full_ajisai()


Open a Terminal: On the offline machine, open a terminal or command prompt. 
Install TensorFlow: Use pip to install the TensorFlow .whl file, specifying the --no-index and --find-links flags. 
pip install --no-index --find-links <path_to_downloaded_files> tensorflow-version.whl 
<path_to_downloaded_files> is the path to the directory where you have the downloaded packages. 
tensorflow-version.whl is the name of the downloaded TensorFlow .whl file. 
Install Dependencies (if needed): If you used pip download, you might have downloaded dependencies. You can install them individually using pip install <package_name.whl>. 
Verify Installation: Check that TensorFlow is installed correctly using pip show tensorflow. 
Important Considerations



https://ieee-dataport.org//open-access/dataset-simulated-inverse-synthetic-aperture-radar-isar-images-automotive-targets
https://essrg.iiitd.edu.in/?page_id=4355
