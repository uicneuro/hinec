#!/usr/bin/env python3
"""
Tractography Slice Viewer GUI

A Python GUI that provides sliders and controls for viewing tractography slices.
Calls the MATLAB visualizeTractographySlices function via subprocess.

Requirements:
    - tkinter (usually included with Python)
    - MATLAB installed and accessible via command line

Usage:
    python tractography_slice_gui.py
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import subprocess
import os
import threading
from pathlib import Path


class TractographySliceGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Tractography Slice Viewer")
        self.root.geometry("600x500")

        # Default values
        self.tracks_file = tk.StringVar(value="tractography_results/tracks_standard.mat")
        self.nim_file = tk.StringVar(value="sample_parcellated.mat")
        self.save_location = tk.StringVar(value="")
        self.x_slice = tk.IntVar(value=64)
        self.y_slice = tk.IntVar(value=64)
        self.z_slice = tk.IntVar(value=32)
        self.tolerance = tk.IntVar(value=2)
        self.show_crosshairs = tk.BooleanVar(value=True)
        self.show_anatomy = tk.BooleanVar(value=True)
        self.color_mode = tk.StringVar(value="direction")
        self.alpha = tk.DoubleVar(value=0.6)

        # Volume dimensions (will be updated when files are loaded)
        self.max_x = 128
        self.max_y = 128
        self.max_z = 64

        self.setup_ui()

    def setup_ui(self):
        """Create the user interface."""
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # File selection section
        file_frame = ttk.LabelFrame(main_frame, text="Data Files", padding="5")
        file_frame.grid(row=0, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))

        ttk.Label(file_frame, text="Tracks File:").grid(row=0, column=0, sticky=tk.W, pady=2)
        ttk.Entry(file_frame, textvariable=self.tracks_file, width=50).grid(row=0, column=1, padx=(5, 5), pady=2)
        ttk.Button(file_frame, text="Browse", command=self.browse_tracks_file).grid(row=0, column=2, pady=2)

        ttk.Label(file_frame, text="NIM File:").grid(row=1, column=0, sticky=tk.W, pady=2)
        ttk.Entry(file_frame, textvariable=self.nim_file, width=50).grid(row=1, column=1, padx=(5, 5), pady=2)
        ttk.Button(file_frame, text="Browse", command=self.browse_nim_file).grid(row=1, column=2, pady=2)

        # Slice controls section
        slice_frame = ttk.LabelFrame(main_frame, text="Slice Positions", padding="5")
        slice_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))

        # X slice (Sagittal)
        ttk.Label(slice_frame, text="X (Sagittal):").grid(row=0, column=0, sticky=tk.W, pady=2)
        x_scale = ttk.Scale(slice_frame, from_=1, to=self.max_x, variable=self.x_slice, orient=tk.HORIZONTAL)
        x_scale.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(5, 5), pady=2)
        x_label = ttk.Label(slice_frame, textvariable=self.x_slice, width=6)
        x_label.grid(row=0, column=2, pady=2)

        # Y slice (Coronal)
        ttk.Label(slice_frame, text="Y (Coronal):").grid(row=1, column=0, sticky=tk.W, pady=2)
        y_scale = ttk.Scale(slice_frame, from_=1, to=self.max_y, variable=self.y_slice, orient=tk.HORIZONTAL)
        y_scale.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=(5, 5), pady=2)
        y_label = ttk.Label(slice_frame, textvariable=self.y_slice, width=6)
        y_label.grid(row=1, column=2, pady=2)

        # Z slice (Axial)
        ttk.Label(slice_frame, text="Z (Axial):").grid(row=2, column=0, sticky=tk.W, pady=2)
        z_scale = ttk.Scale(slice_frame, from_=1, to=self.max_z, variable=self.z_slice, orient=tk.HORIZONTAL)
        z_scale.grid(row=2, column=1, sticky=(tk.W, tk.E), padx=(5, 5), pady=2)
        z_label = ttk.Label(slice_frame, textvariable=self.z_slice, width=6)
        z_label.grid(row=2, column=2, pady=2)

        # Store scale references for updating ranges
        self.x_scale = x_scale
        self.y_scale = y_scale
        self.z_scale = z_scale

        # Configure column weights
        slice_frame.columnconfigure(1, weight=1)

        # Options section
        options_frame = ttk.LabelFrame(main_frame, text="Display Options", padding="5")
        options_frame.grid(row=2, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))

        # Tolerance
        ttk.Label(options_frame, text="Slice Thickness:").grid(row=0, column=0, sticky=tk.W, pady=2)
        ttk.Spinbox(options_frame, from_=1, to=10, textvariable=self.tolerance, width=10).grid(row=0, column=1, sticky=tk.W, padx=(5, 0), pady=2)
        ttk.Label(options_frame, text="voxels").grid(row=0, column=2, sticky=tk.W, padx=(5, 0), pady=2)

        # Color mode
        ttk.Label(options_frame, text="Color Mode:").grid(row=1, column=0, sticky=tk.W, pady=2)
        color_combo = ttk.Combobox(options_frame, textvariable=self.color_mode, values=["direction", "uniform"], width=12)
        color_combo.grid(row=1, column=1, sticky=tk.W, padx=(5, 0), pady=2)
        color_combo.state(['readonly'])

        # Alpha (transparency)
        ttk.Label(options_frame, text="Background Alpha:").grid(row=2, column=0, sticky=tk.W, pady=2)
        alpha_scale = ttk.Scale(options_frame, from_=0.0, to=1.0, variable=self.alpha, orient=tk.HORIZONTAL, length=150)
        alpha_scale.grid(row=2, column=1, sticky=tk.W, padx=(5, 0), pady=2)
        alpha_label = ttk.Label(options_frame, text=f"{self.alpha.get():.1f}", width=6)
        alpha_label.grid(row=2, column=2, padx=(5, 0), pady=2)

        # Update alpha label when value changes
        def update_alpha_label(*args):
            alpha_label.config(text=f"{self.alpha.get():.1f}")
        self.alpha.trace('w', update_alpha_label)

        # Checkboxes
        ttk.Checkbutton(options_frame, text="Show crosshairs", variable=self.show_crosshairs).grid(row=3, column=0, columnspan=2, sticky=tk.W, pady=2)
        ttk.Checkbutton(options_frame, text="Show anatomy background", variable=self.show_anatomy).grid(row=4, column=0, columnspan=2, sticky=tk.W, pady=2)

        # Save section
        save_frame = ttk.LabelFrame(main_frame, text="Output", padding="5")
        save_frame.grid(row=3, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))

        ttk.Label(save_frame, text="Save Location:").grid(row=0, column=0, sticky=tk.W, pady=2)
        ttk.Entry(save_frame, textvariable=self.save_location, width=50).grid(row=0, column=1, padx=(5, 5), pady=2)
        ttk.Button(save_frame, text="Browse", command=self.browse_save_location).grid(row=0, column=2, pady=2)

        save_frame.columnconfigure(1, weight=1)

        # Action buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=4, column=0, columnspan=2, pady=(10, 0))

        ttk.Button(button_frame, text="Preview", command=self.preview_slices, style="Accent.TButton").grid(row=0, column=0, padx=(0, 10))
        ttk.Button(button_frame, text="Render & Save", command=self.render_and_save).grid(row=0, column=1, padx=(0, 10))
        ttk.Button(button_frame, text="Reset to Center", command=self.reset_to_center).grid(row=0, column=2)

        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
        status_bar.grid(row=5, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(10, 0))

        # Configure main grid weights
        main_frame.columnconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

    def browse_tracks_file(self):
        """Browse for tracks file."""
        filename = filedialog.askopenfilename(
            title="Select Tracks File",
            filetypes=[("MAT files", "*.mat"), ("All files", "*.*")]
        )
        if filename:
            self.tracks_file.set(filename)

    def browse_nim_file(self):
        """Browse for NIM file."""
        filename = filedialog.askopenfilename(
            title="Select NIM File",
            filetypes=[("MAT files", "*.mat"), ("All files", "*.*")]
        )
        if filename:
            self.nim_file.set(filename)

    def browse_save_location(self):
        """Browse for save location."""
        filename = filedialog.asksaveasfilename(
            title="Save Slice View As",
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")]
        )
        if filename:
            self.save_location.set(filename)

    def reset_to_center(self):
        """Reset slice positions to center of volume."""
        self.x_slice.set(self.max_x // 2)
        self.y_slice.set(self.max_y // 2)
        self.z_slice.set(self.max_z // 2)
        self.status_var.set("Reset to center positions")

    def validate_files(self):
        """Validate that required files exist."""
        tracks_path = self.tracks_file.get()
        nim_path = self.nim_file.get()

        if not tracks_path or not os.path.exists(tracks_path):
            messagebox.showerror("Error", f"Tracks file not found: {tracks_path}")
            return False

        if not nim_path or not os.path.exists(nim_path):
            messagebox.showerror("Error", f"NIM file not found: {nim_path}")
            return False

        return True

    def build_matlab_command(self, save_output=False):
        """Build the MATLAB command to execute."""
        tracks_file = self.tracks_file.get()
        nim_file = self.nim_file.get()
        x = self.x_slice.get()
        y = self.y_slice.get()
        z = self.z_slice.get()

        # Build base command
        cmd_parts = [
            f"visualizeTractographySlices('{tracks_file}', '{nim_file}', {x}, {y}, {z}",
            f"'tolerance', {self.tolerance.get()}",
            f"'show_crosshairs', {str(self.show_crosshairs.get()).lower()}",
            f"'show_anatomy', {str(self.show_anatomy.get()).lower()}",
            f"'color_mode', '{self.color_mode.get()}'",
            f"'alpha', {self.alpha.get()}"
        ]

        # Add save parameter if needed
        if save_output and self.save_location.get():
            cmd_parts.append(f"'save', '{self.save_location.get()}'")

        matlab_cmd = ", ".join(cmd_parts) + ");"

        return matlab_cmd

    def run_matlab_command(self, matlab_cmd, success_msg="Command completed"):
        """Run MATLAB command in a separate thread."""
        def run_in_thread():
            try:
                self.status_var.set("Running MATLAB command...")
                self.root.update()

                # Execute MATLAB command
                full_cmd = ["matlab", "-batch", matlab_cmd]
                result = subprocess.run(full_cmd, capture_output=True, text=True, timeout=300)

                if result.returncode == 0:
                    self.status_var.set(success_msg)
                else:
                    error_msg = f"MATLAB error: {result.stderr}"
                    self.status_var.set("Error occurred")
                    messagebox.showerror("MATLAB Error", error_msg)

            except subprocess.TimeoutExpired:
                self.status_var.set("MATLAB command timed out")
                messagebox.showerror("Error", "MATLAB command timed out after 5 minutes")

            except Exception as e:
                self.status_var.set("Error occurred")
                messagebox.showerror("Error", f"Failed to run MATLAB command: {str(e)}")

        # Run in separate thread to prevent GUI freezing
        thread = threading.Thread(target=run_in_thread)
        thread.daemon = True
        thread.start()

    def preview_slices(self):
        """Preview the slice view (display only)."""
        if not self.validate_files():
            return

        matlab_cmd = self.build_matlab_command(save_output=False)
        self.run_matlab_command(matlab_cmd, "Preview displayed")

    def render_and_save(self):
        """Render and save the slice view."""
        if not self.validate_files():
            return

        if not self.save_location.get():
            messagebox.showerror("Error", "Please specify a save location")
            return

        matlab_cmd = self.build_matlab_command(save_output=True)
        save_path = self.save_location.get()
        self.run_matlab_command(matlab_cmd, f"Saved to {os.path.basename(save_path)}")


def main():
    """Main entry point."""
    root = tk.Tk()
    app = TractographySliceGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()