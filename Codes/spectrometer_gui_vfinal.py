# Python application to interface with an Arduino board for capturing spectra
# using the HAMAMATSU C12880MA mini-spectrometer.
#
# Designed for use on a Raspberry Pi.
#
# This code was developed with the invaluable assistance of DeepSeek (https://chat.deepseek.com/),
# which provided not only technical expertise but also advanced Spanish language support,
# enabling clear collaboration even for complex concepts.

# Import required libraries
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import serial
import numpy as np
from scipy.io import savemat
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import Image, ImageTk
import os
import time

# Global Constants
# ++++++++++++++++++
# Serial port configuration
SERIAL_BAUDRATE = 115200  # Baud rate for Arduino communication

# Spectrometer configuration
DATA_POINTS = 288  # Number of spectral channels from C12880MA

# Calibration coefficients (λ(nm) = C0 + C1*n + C2*n² + ...)
A0 = 3.213544722E+02
B1 = 2.687903253E+00
B2 = -9.427462678E-04
B3 = -9.880332677E-06
B4 = 1.833932084E-08
B5 = -8.146778297E-12

# Calibration vector
CALIBRATION = [A0, B1, B2, B3, B4, B5]

# UI configuration
font_size = 12  # Base font size for UI elements

class SpectrometerApp:
    def __init__(self, root):
        """Initialize the spectrometer application.
        
        Args:
            root (tk.Tk): The main Tkinter window.
        """
        # Store reference to root window
        self.root = root
        
        # Configure main window
        self.root.title("C12880MA Spectrometer Controller")
        self.root.update()  # Force initial UI update

        # State variables
        self.vlines = []  # Stores vertical reference lines on plots
        self.state_oscilloscope = 0  # Tracks oscilloscope mode state
        self.serial_port = None  # Serial port object
        
        # Data storage arrays
        self.dark_data = np.zeros(DATA_POINTS)  # Dark reference spectrum
        self.ref_data = np.zeros(DATA_POINTS)    # Reference spectrum
        self.raw_data = np.zeros(DATA_POINTS)    # Raw measured spectrum
        self.continuous_buffer = np.zeros(DATA_POINTS)  # Buffer for averaging

        # Measurement parameters
        self.measure_count = 0  # Spectrum counter
        self.delay = 1         # Integration time in microseconds
        self.iterations = 10   # Default number of measurement iterations

        # Wavelength calibration
        # Note: The -20 offset was determined empirically during lab calibration
        self.wavelengths = np.polyval(CALIBRATION[::-1], np.arange(DATA_POINTS)) - 20 
        
        # Plotting flags
        self.plotting_avg_min = 1  # Toggles minimum wavelength plotting

        # Configure UI grid weights
        self.root.grid_columnconfigure(1, weight=1)  # Expand graph column
        self.root.grid_rowconfigure(0, weight=1)     # Expand logo row
        self.root.grid_rowconfigure(1, weight=0)     # Fixed height for controls
        self.root.grid_rowconfigure(2, weight=0)
        self.root.grid_rowconfigure(3, weight=0)
        self.root.grid_rowconfigure(4, weight=0)

        # Initialize UI and serial connection
        self.setup_ui()
        self.init_serial()

    def setup_ui(self):
        """Configure all user interface elements."""
        
        # --- Operations Frame (row 1, column 0) ---
        ops_frame = tk.LabelFrame(self.root, text="Operations", font=("Arial", font_size))
        ops_frame.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")
        ops_frame.grid_columnconfigure(0, weight=1)  # Expand column
        
        # Inner frame for consistent background
        frame_interior_ops_frame = tk.Frame(
            ops_frame, 
            bg="dark sea green",
            bd=0,
            highlightthickness=0  
        )
        frame_interior_ops_frame.grid(row=0, column=0, sticky="nsew", padx=0, pady=0)
        frame_interior_ops_frame.grid_columnconfigure(0, weight=1)

        # Control buttons
        self.btn_dark = tk.Button(frame_interior_ops_frame, text="DARK", 
                                command=self.capture_dark, state="disabled", 
                                font=("Arial", font_size))
        self.btn_dark.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        self.btn_ref = tk.Button(frame_interior_ops_frame, text="REF.", 
                               command=self.capture_reference, state="disabled", 
                               font=("Arial", font_size))
        self.btn_ref.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")

        self.btn_data = tk.Button(frame_interior_ops_frame, text="DATA", 
                                command=self.capture_data, state="disabled", 
                                font=("Arial", font_size))
        self.btn_data.grid(row=2, column=0, padx=5, pady=5, sticky="nsew")
        
        self.btn_oscill = tk.Button(frame_interior_ops_frame, text="OSCILLOSCOPE", 
                                  command=self.mode_oscilloscope, state="disabled", 
                                  font=("Arial", font_size))
        self.btn_oscill.grid(row=3, column=0, padx=5, pady=5, sticky="nsew")
        
        # --- Settings Frame (row 2, column 0) ---
        settings_frame = tk.LabelFrame(self.root, text="Settings", font=("Arial", font_size))
        settings_frame.grid(row=2, column=0, padx=5, pady=5, sticky="nsew")
        settings_frame.grid_columnconfigure(0, weight=1)
        settings_frame.grid_rowconfigure(0, weight=1)

        frame_interior_settings_frame = tk.Frame(
            settings_frame, 
            bg="dark sea green",
            bd=0,
            highlightthickness=0
        )
        frame_interior_settings_frame.grid(row=0, column=0, sticky="nsew", padx=0, pady=0)
        frame_interior_settings_frame.grid_columnconfigure(0, weight=1)

        # Iterations control
        label_iterations = tk.Label(frame_interior_settings_frame, text="Iterations", 
                                  font=("Arial", font_size), bg="dark sea green")
        label_iterations.grid(row=0, column=0, padx=5, pady=5, sticky="w")

        self.entry_iter = tk.Entry(frame_interior_settings_frame, font=("Arial", font_size), width=5)
        self.entry_iter.insert(0, "15")
        self.entry_iter.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")

        # Integration time control
        label_delay = tk.Label(frame_interior_settings_frame, text="Integ. time (μs):", 
                              font=("Arial", font_size), bg="dark sea green")
        label_delay.grid(row=2, column=0, padx=5, pady=5, sticky="w")

        self.entry_delay = tk.Entry(frame_interior_settings_frame, font=("Arial", font_size), width=5)
        self.entry_delay.insert(0, "1")
        self.entry_delay.grid(row=3, column=0, padx=5, pady=5, sticky="nsew")

        self.btn_setting = tk.Button(frame_interior_settings_frame, text="SET", 
                                   command=self.apply_settings, font=("Arial", font_size))
        self.btn_setting.grid(row=4, column=0, padx=5, pady=5, sticky="nsew")
               
        # --- Find Minimum Frame (row 3, column 0) ---
        findmin_frame = tk.LabelFrame(self.root, text="Find Minimum", font=("Arial", font_size))
        findmin_frame.grid(row=3, column=0, padx=5, pady=5, sticky="nsew")
        findmin_frame.grid_columnconfigure(0, weight=1)

        frame_interior_findmin_frame = tk.Frame(
            findmin_frame, 
            bg="dark sea green",
            bd=0,
            highlightthickness=0
        )
        frame_interior_findmin_frame.grid(row=0, column=0, sticky="nsew", padx=0, pady=0)
        frame_interior_findmin_frame.grid_columnconfigure(0, weight=1)

        # Wavelength range controls
        label_wl_lower = tk.Label(frame_interior_findmin_frame, text="Wl_lower (nm)", 
                                font=("Arial", font_size), bg="dark sea green")
        label_wl_lower.grid(row=0, column=0, padx=5, pady=5, sticky="w")

        self.entry_wl_lower = tk.Entry(frame_interior_findmin_frame, font=("Arial", font_size), width=5)
        self.entry_wl_lower.insert(0, "450.0")
        self.entry_wl_lower.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")

        label_wl_upper = tk.Label(frame_interior_findmin_frame, text="Wl_upper (nm)", 
                                font=("Arial", font_size), bg="dark sea green")
        label_wl_upper.grid(row=2, column=0, padx=5, pady=5, sticky="w")

        self.entry_wl_upper = tk.Entry(frame_interior_findmin_frame, font=("Arial", font_size), width=5)
        self.entry_wl_upper.insert(0, "650.0")
        self.entry_wl_upper.grid(row=3, column=0, padx=5, pady=5, sticky="nsew")

        # Minimum finder button
        self.btn_findmin = tk.Button(frame_interior_findmin_frame, width=5, text="FIND MIN.", 
                                   command=self.find_minimum, state="disabled", 
                                   font=("Arial", font_size))
        self.btn_findmin.grid(row=4, column=0, padx=5, pady=5, sticky="nsew")

        # --- Progress Bar (row 4, columns 0-1) ---
        self.progress_label = tk.Label(self.root, text="Progress:")
        self.progress_label.grid(row=4, column=0, padx=5, pady=(0, 0), sticky="e")

        self.progress = ttk.Progressbar(self.root, orient="horizontal", mode="determinate")
        self.progress.grid(row=4, column=1, padx=0, pady=(0, 0), sticky="ew")

        # --- Graph Frame (column 1, rows 0-3) ---
        graph_frame = ttk.Frame(self.root)
        graph_frame.grid(row=0, column=1, rowspan=4, padx=0, pady=0, sticky="nsew")

        # Configure matplotlib figure
        self.fig = Figure(figsize=(6, 4), dpi=80)

        # Top subplot: Spectra visualization
        self.ax = self.fig.add_subplot(211)
        self.ax.set_ylim(0.0, 5.0)
        self.ax.set_xlim(self.wavelengths[0], self.wavelengths[-1])
        
        # Bottom subplot: Minimum wavelength tracking
        self.ax1 = self.fig.add_subplot(212)
        
        # Configure plot lines
        self.line_avg_min, = self.ax1.plot([], [],
           color=(0, 0, 1, 0.25),  # Semi-transparent blue line
           marker='o',             # Circular markers
           markerfacecolor='red',  # Red fill
           visible=False
        )
        
        # Set initial time axis limits (2.4s per spectrum)
        self.ax1.set_xlim(0.0, 2.4*5)
                
        # Normalized spectrum line
        self.line_Normalized,  = self.ax.plot([], [], 'g-', linewidth=1.5, 
                                             label="Normalized data", visible=False)
        
        # Minimum point marker
        self.point_min, = self.ax.plot([], [], 'ro', markersize=8, 
                                      label='Minimum', visible=False)

        # Minimum wavelength label
        self.text_minimum = self.ax.text(
                0, 0,
                "",
                verticalalignment='bottom',
                fontsize=8,
                color='red'
            )
        self.text_minimum.set_visible(False)
        
        # Configure auto-scaling
        self.ax.set_autoscale_on(True)
        self.ax.autoscale(enable=True, axis='y')
        self.ax1.set_autoscale_on(True)
        self.ax1.autoscale(enable=True, axis='y')
        
        # Configure grid lines
        self.ax.grid(True, axis='x', color='lightgray', linestyle='--', linewidth=0.8)
        self.ax.grid(True, axis='y', color='lightgray', linestyle='--', linewidth=0.8)
        self.ax1.grid(True, axis='x', color='lightgray', linestyle='--', linewidth=0.8)
        self.ax1.grid(True, axis='y', color='lightgray', linestyle='--', linewidth=0.8)
        
        # Axis labels
        self.ax.set_xlabel("Wavelength (nm)", fontsize=12)
        self.ax.set_ylabel("Intensity (V)", fontsize=12)
        self.ax1.set_ylabel("Wavelength (nm)", fontsize=12)
        self.ax1.set_xlabel("Minimum finder: Elapsed time (s)", fontsize=12)
        self.ax.tick_params(axis='both', labelsize=12)

        # Set fixed wavelength ticks
        self.ax.set_xticks([325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 
                          575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 
                          825, 850, 875, 900])

        # Initialize plot lines
        self.line_Average,  = self.ax.plot([], [], 'b-', linewidth=1.5, 
                                         label="Avg. data", visible=False)
        self.line_Raw_data, = self.ax.plot([], [], 'r-', alpha=0.3, 
                                         label="Raw data", visible=False)
        self.point_min_avg, = self.ax.plot([], [], 'ro', markersize=8, 
                                         visible=False)
        
        # Adjust figure margins
        self.fig.subplots_adjust(
            left=0.06,     # 6% left margin
            right=0.95,    # 5% right margin
            bottom=0.07,   # 7% bottom margin
            top=0.97       # 3% top margin
        )

        # Embed matplotlib figure in Tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=graph_frame)
        self.canvas.get_tk_widget().pack(expand=True, fill="both")

        # --- Application Logo (row 0, column 0) ---
        self.logo_frame = tk.Frame(self.root)
        self.logo_frame.grid(row=0, column=0, sticky="nsew", padx=1, pady=1)
        self.logo_frame.config(relief="solid", borderwidth=1)

        # Load and resize logo image
        try:
            image = Image.open("logo_UB.jpg")
            original_width, original_height = image.size
            
            # Force UI update to obtain correct values for frame_height and frame_width
            self.root.update()
            
            # Calculate dimensions maintaining aspect ratio
            frame_height = self.logo_frame.winfo_height()
            aspect_ratio = original_height / original_width
            frame_width = int(frame_height / aspect_ratio)
            
            # Resize using best available filter
            try:
                resize_filter = Image.Resampling.LANCZOS  # Pillow >= 9.1.0
            except AttributeError:
                resize_filter = Image.LANCZOS  # Older versions

            resized_image = image.resize((frame_width, frame_height), resize_filter)
            self.logo_img = ImageTk.PhotoImage(resized_image)

            self.logo_label = tk.Label(self.logo_frame, image=self.logo_img, bg='white')
            self.logo_label.pack(expand=True, fill="both")
        except Exception as e:
            print(f"Error loading logo: {e}")

        # Show keyboard shortcuts info
        messagebox.showinfo("Information", 
                          "Keyboard Shortcuts:\n"
                          "- (m/M): Save last spectrum to TXT\n"
                          "- (s/S): Toggle minimum wavelength plotting")
        
        # Bind keyboard events
        self.root.bind("m", self.save_TXT)
        self.root.bind("M", self.save_TXT)
        self.root.bind("s", self.enable_plotting_avg_min)
        self.root.bind("S", self.enable_plotting_avg_min)
        
        # Initialize save TXT flag
        self.enable_save_TXT = False
        
        # Force UI update
        self.root.update()
    
    def enable_plotting_avg_min(self, event):
        """Toggle minimum wavelength plotting on/off.
        
        Triggered by 's' or 'S' key press.
        """
        self.plotting_avg_min ^= 1  # Toggle between 0 and 1

    def save_TXT(self, event):
        """Save raw spectrum data to text file.
        
        Triggered by 'm' or 'M' key press when enabled.
        """
        if self.enable_save_TXT:
            self.enable_save_TXT = False  # Disable after saving
            
            # Create 2D array of wavelengths and intensities
            data = np.column_stack((self.wavelengths, self.raw_data))

            # Prompt for save location
            filename = filedialog.asksaveasfilename(defaultextension=".txt")
            if filename:
                np.savetxt(filename, data, delimiter=",")            
        
    def update_legend(self):
        """Update the plot legend with currently visible lines."""
        self.canvas.draw()
        self.root.update()
        
        # Collect visible lines and labels from both axes
        lines = []
        labels = []
    
        # Axis lines
        for line in self.ax.get_lines():
            if line.get_visible():
                label = line.get_label()
                if not label.startswith('_'):
                    lines.append(line)
                    labels.append(label)
        
        # Remove existing legend if present
        if hasattr(self, 'legend'):
            self.legend.remove()
    
        # Create new unified legend
        self.legend = self.ax.legend(
            lines, 
            labels,
            loc='best'  # Automatic positioning
        )
        
        # Update legend in figure
        self.canvas.draw_idle()
        self.root.update()

    def mode_oscilloscope(self):
        """Toggle oscilloscope mode for hardware debugging."""
        if not self.state_oscilloscope:
            # Enter oscilloscope mode
            self.clear_plots()
            
            # Send command to Arduino
            self.serial_port.write(b"OSCILLOSCOPE\n")
            
            # Update UI state
            self.state_oscilloscope = 1
            self.btn_oscill.config(text="STOP OSCILL.")
            
            # Disable other controls
            self.toggle_buttons_state("disabled")
            self.canvas.draw()
            self.root.update()

        else:
            # Exit oscilloscope mode
            self.serial_port.write(b"STOP OSCILLOSCOPE\n")
            
            # Wait for confirmation from Arduino
            while True:
                if self.serial_port.in_waiting > 0:
                    data = self.serial_port.readline().decode().strip()
                    if data == "STOP OSCILLOSCOPE":
                        break
                    else:
                        messagebox.showerror("ERROR", "Failed to exit oscilloscope mode")
                        self.root.destroy()
                else:
                    time.sleep(0.1)
        
            # Restore UI state
            self.state_oscilloscope = 0
            self.btn_oscill.config(text="OSCILLOSCOPE")
            self.restore_button_states()
            self.canvas.draw()
            self.root.update()
 
    def remove_vlines(self):
        """Remove all vertical reference lines from plots."""
        for line in self.vlines:
            line.remove()
        self.vlines.clear()
                   
    def find_minimum(self):
        """Find and mark the minimum intensity in the specified wavelength range."""
        try:
            # Reset progress bar
            self.progress['value'] = 0
            self.progress_label.config(text="Progress:")
            self.root.update()

            # Clear existing plots
            self.clear_plots()
            self.line_Normalized.set_visible(True)
            self.point_min.set_visible(True)            
            self.remove_vlines()
        
            # Validate wavelength range
            wl_lower = float(self.entry_wl_lower.get())
            wl_upper = float(self.entry_wl_upper.get())

            if wl_lower >= wl_upper:
                raise ValueError("Lower limit must be less than upper limit")
            if wl_lower < self.wavelengths[0] or wl_upper > self.wavelengths[-1]:
                raise ValueError(f"Valid range: {self.wavelengths[0]:.1f} nm to {self.wavelengths[-1]:.1f} nm")

            # Plot normalized data
            self.line_Normalized.set_data(self.wavelengths, self.processed_data)
            self.ax.relim()
            self.ax.autoscale_view(scaley=True)

            # Add range boundary lines
            self.vlines.append(self.ax.axvline(x=wl_lower, color='black', alpha=0.95, linewidth=2))
            self.vlines.append(self.ax.axvline(x=wl_upper, color='black', alpha=0.95, linewidth=2))

            # Find minimum in specified range
            mask = (self.wavelengths >= wl_lower) & (self.wavelengths <= wl_upper)
            if not np.any(mask):
                raise ValueError("No data in specified range")

            min_idx = np.argmin(self.processed_data[mask])
            min_wl = self.wavelengths[mask][min_idx]
            min_intensity = self.processed_data[mask][min_idx]

            # Mark minimum with red dot and label
            self.point_min.set_data(min_wl, min_intensity)
            self.text_minimum = self.ax.text(
                min_wl, min_intensity-0.05,
                f'  Minimum: {min_wl:.1f} nm\n  Intensity: {min_intensity:.4f}',
                verticalalignment='bottom',
                fontsize=8,
                color='red'
            )

            # Update ploting
            self.update_legend()
            self.canvas.draw()

        except ValueError as e:
            self.ax.set_title(f"Error: {str(e)}", fontsize=10, color='red')
            self.canvas.draw()
            print(f"Error: {e}")

    def init_serial(self):
        """Initialize serial connection to Arduino."""
        # Try common Linux serial ports
        port_device = "/dev/ttyUSB0" if os.path.exists("/dev/ttyUSB0") else "/dev/ttyUSB1"
        
        try:
            self.serial_port = serial.Serial(port_device, SERIAL_BAUDRATE, timeout=2)
            print(f"Connected to Arduino at {port_device}")
        except Exception as e:
            print(f"Serial error: {e}")
            messagebox.showerror("Serial Error", f"Failed to connect: {e}")

    def apply_settings(self):
        """Apply measurement settings from UI to Arduino."""
        try:
            # Get values from UI
            self.iterations = int(self.entry_iter.get())
            self.ax1.set_xlim(0.0, 2.4*self.iterations)  # Adjust time axis
            
            # Only allow delay setting once per session
            if self.entry_delay["state"] == "normal":
                self.delay = int(self.entry_delay.get())

            # Configure progress bar
            self.progress['maximum'] = self.iterations
            self.progress['value'] = 0

            # Send settings to Arduino
            if self.serial_port:
                command = f"SETTING:{self.delay},{self.iterations}\n"
                expected_response = f"READY:{self.delay},{self.iterations}"
                
                self.serial_port.write(command.encode())
                response = self.serial_port.readline().decode().strip()
                print(response)

                if response == expected_response:
                    self.btn_dark["state"] = "normal"
                    self.btn_oscill["state"] = "normal"
                    self.entry_delay.configure(state="disabled", 
                                             disabledbackground="light gray")
                    print("Settings applied successfully")
                else:
                    print("Error: Invalid Arduino response")
        
        except ValueError:
            messagebox.showerror("ERROR", "Invalid values in Delay/Iterations")
        
        except Exception as e:
            messagebox.showerror("ERROR", f"Settings error: {e}")

    def capture_spectrum(self, mode):
        """Capture and process spectra based on measurement mode.
        
        Args:
            mode (str): Measurement mode ("DARK", "REFERENCE", or "DATA")
        """
        # Disable UI during capture
        self.toggle_buttons_state("disabled")
        
        # Delete previous vertical lines and text for minimum
        self.remove_vlines()
        self.text_minimum.set_visible(False)
        
        # Initialize tracking arrays
        self.x_min_avg = [None]*self.iterations
        self.y_min_avg = [None]*self.iterations
        self.line_avg_min.set_data([], [])
        self.line_avg_min.set_visible(True)

        try:
            # Mode-specific configuration
            title_map = {
                    "DARK": "DARK Acquisition",
                    "REFERENCE": "REFERENCE Acquisition",
                    "DATA": "DATA Acquisition"
            }

            # Reset progress
            self.progress['value'] = 0

            if not self.serial_port:
                return

            # Clear existing plots
            self.clear_plots()
            
            # Configure lines 
            self.line_Average.set_visible(True)
            self.point_min_avg.set_visible(True)
            
            if mode in ["DARK", "REFERENCE"]:
                self.line_Raw_data.set_visible(True)
                self.ax.set_ylabel("Wavelength (nm)", fontsize=12)
        
            else:
                self.line_Normalized.set_visible(True)
                self.ax.set_ylabel("Intensity normalized", fontsize=12)                
            
            # Title for the plot
            self.ax.set_title(title_map[mode], fontsize=12)
            
            # Initialize data structures
            self.continuous_buffer = np.zeros(DATA_POINTS)
            self.measure_count = 0
            self.progress["value"] = 0
            self.progress_label.config(text="Progress: 0")

            # Update plotting area
            self.update_legend()
            self.canvas.draw_idle()
            self.root.update()
            
            # Set time for the capture start            
            start_time = time.time()
            
            # Send command to Arduino: Start acquisition
            self.serial_port.write(b"GET DATA\n")

            while True:
                line = self.serial_port.readline().decode().strip()

                if line.startswith("DATA:"):
                    # Get elapsed time since start
                    elapsed = time.time() - start_time
                    
                    # Increment measuring counter
                    self.measure_count += 1
          
                    # Convert raw ADC values to voltage (0-5V)
                    self.raw_data = 5.0 * np.array(list(map(int, line.split(":")[1].split(",")))) / 1023.0
                
                    # Update data plot
                    if mode in ["DARK", "REFERENCE"]:   
                        # Accumulate for averaging
                        self.continuous_buffer += self.raw_data
                             
                        # Update raw data plot
                        self.line_Raw_data.set_data(self.wavelengths, self.raw_data)
                        
                    else:
                        self.processed_data = (self.raw_data - self.dark_data) / self.ref_data
                        
                        # Accumulate for averaging
                        self.continuous_buffer += self.processed_data
                        
                        # Update normalized data plot
                        self.line_Normalized.set_data(self.wavelengths, self.processed_data)
                        
                    # Calculate running average
                    avg_data = self.continuous_buffer / self.measure_count                    
                    
                    # Update running average plot
                    self.line_Average.set_data(self.wavelengths, avg_data)    
                    
                    # Track minimum wavelength
                    min_idx = np.argmin(avg_data)
                    min_wl = self.wavelengths[min_idx]
                    self.point_min_avg.set_data(min_wl, avg_data[min_idx])
                             
                    self.x_min_avg[self.measure_count-1] = elapsed
                    self.y_min_avg[self.measure_count-1] = min_wl
                    
                    # Update minimum tracking plot if enabled
                    if self.plotting_avg_min:
                        self.line_avg_min.set_data(self.x_min_avg, self.y_min_avg)
                        self.ax1.relim()
                        self.ax1.autoscale_view(scaley=True)
                                            
                    # Update progress
                    self.progress["value"] += 1
                    self.progress_label.config(text=f"Progress: {self.progress['value']}")
                    
                    # Enable data saving
                    self.enable_save_TXT = True 
                    
                    # Refresh plots
                    self.ax.set_title(title_map[mode], fontsize=12)
                    self.ax.relim()
                    self.ax.autoscale_view(scaley=True)
                    
                    self.update_legend()
                    self.canvas.draw()
                    self.root.update()
                    
                    # Request next spectrum to Arduino
                    time.sleep(0.00001)  # Brief pause
                    self.serial_port.write(b"NEXT\n")                    

                elif line == "END DATA":
                    self.progress_label.config(text="Finished!!!")

                    # Save minimum tracking data
                    data = np.column_stack((np.array(self.x_min_avg), np.array(self.y_min_avg)))
                    filename = f"minimum_averaged_{time.strftime('%Y_%m_%d-%H_%M_%S')}.txt"
                    np.savetxt(filename, data, delimiter=",") 
                
                    # Final processing
                    averaged_data = self.continuous_buffer / self.measure_count
                    self.continuous_buffer = np.zeros(DATA_POINTS)
                    self.measure_count = 0

                    # Mode-specific processing
                    if mode == "DARK":
                        self.dark_data = averaged_data
                        
                    elif mode == "REFERENCE":
                        self.ref_data = averaged_data - self.dark_data
                        
                        # To avoid division by zero errors, we substitute 0.0001 for the possible zeros in the matrix self.ref_data
                        self.ref_data = np.where(self.ref_data == 0, 0.001, self.ref_data)
                        
                    elif mode == "DATA":
                        self.processed_data = (averaged_data - self.dark_data) / self.ref_data
                   
                        # Save processed data
                        filename = filedialog.asksaveasfilename(defaultextension=".mat")
                        if filename:
                            savemat(filename, {
                                'wavelengths': self.wavelengths,
                                'spectrum': self.processed_data,
                                'dark_spectrum': self.dark_data,
                                'ref_spectrum': self.ref_data
                            })

                    break  # Exit acquisition loop

        except Exception as e:
            print(f"Capture error: {e}")
            self.ax.set_title(f"Error: {str(e)}", fontsize=12, color='red')
            self.canvas.draw()

        finally:
            # Always restore UI state
            self.restore_button_states(mode)

    def capture_dark(self):
        """Capture dark reference spectrum."""
        self.capture_spectrum("DARK")
        self.btn_ref["state"] = "normal"  # Enable reference capture

    def capture_reference(self):
        """Capture reference spectrum."""
        self.capture_spectrum("REFERENCE")
        self.btn_data["state"] = "normal"  # Enable data capture

    def capture_data(self):
        """Capture sample spectrum."""
        self.capture_spectrum("DATA")
        self.btn_findmin["state"] = "normal"  # Enable minimum finding

    def toggle_buttons_state(self, state="normal"):
        """Enable/disable all control buttons."""
        buttons = [
            self.btn_dark,
            self.btn_ref,
            self.btn_data,
            self.btn_findmin,
            self.btn_setting,
            self.btn_oscill
        ]
        for btn in buttons:
            btn.config(state=state)

    def restore_button_states(self, mode=None):
        """Restore button states based on current mode."""
        self.btn_setting.config(state="normal")
        self.btn_dark.config(state="normal")
        self.btn_oscill.config(state="normal")
        
        if mode == "DARK":
            self.btn_ref.config(state="normal")
        elif mode == "REFERENCE":
            self.btn_ref.config(state="normal")
            self.btn_data.config(state="normal")
        elif mode == "DATA":
            self.btn_ref.config(state="normal")
            self.btn_data.config(state="normal")
            self.btn_findmin.config(state="normal")

    def clear_plots(self):
        """Clear all plot data."""
        self.line_Average.set_data([], [])
        self.line_Raw_data.set_data([], [])
        self.line_Normalized.set_data([], [])
        self.point_min.set_data([], [])
        self.point_min_avg.set_data([], [])
        self.line_avg_min.set_data([], [])
        self.line_Average.set_visible(False)
        self.line_Raw_data.set_visible(False)
        self.line_Normalized.set_visible(False)
        self.point_min.set_visible(False)
        self.text_minimum.set_visible(False)
        self.update_legend()

if __name__ == "__main__":
    # Create and configure main window
    root = tk.Tk()
    root.geometry(f"{root.winfo_screenwidth()}x{root.winfo_screenheight()}")
    root.update()
            
    # Start application
    app = SpectrometerApp(root)
    root.mainloop()
