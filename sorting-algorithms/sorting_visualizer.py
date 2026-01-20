# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sounddevice as sd

# Configuration
NUM_BARS = 20
DELAY_MS = 80  # Slightly increased delay to allow sound to play
SOUND_ENABLED = True
SAMPLE_RATE = 44100

# Sound parameters
MIN_FREQ = 200   # Hz for smallest difference
MAX_FREQ = 1200  # Hz for largest difference
BEEP_DURATION = 0.04  # seconds

# Generate random data
np.random.seed(42)
original_data = np.random.randint(1, 100, NUM_BARS)
max_diff = 100  # Maximum possible difference

def generate_beep(frequency, duration=BEEP_DURATION):
    """Generate a beep waveform"""
    t = np.linspace(0, duration, int(SAMPLE_RATE * duration), False)
    # Generate sine wave with envelope (fade in/out to avoid clicks)
    envelope = np.ones_like(t)
    fade_samples = int(SAMPLE_RATE * 0.005)  # 5ms fade
    if fade_samples > 0:
        envelope[:fade_samples] = np.linspace(0, 1, fade_samples)
        envelope[-fade_samples:] = np.linspace(1, 0, fade_samples)
    wave = 0.3 * np.sin(2 * np.pi * frequency * t) * envelope
    return wave.astype(np.float32)

def play_beep(frequency):
    """Play a beep sound - stops previous sound first"""
    if not SOUND_ENABLED:
        return
    try:
        sd.stop()  # Stop any currently playing sound
        wave = generate_beep(frequency)
        sd.play(wave, SAMPLE_RATE)
    except Exception:
        pass  # Silently ignore audio errors

def get_frequency_for_diff(diff):
    """Map value difference to frequency"""
    normalized = min(abs(diff) / max_diff, 1.0)
    return MIN_FREQ + normalized * (MAX_FREQ - MIN_FREQ)

# Store all operations for each sorting algorithm
def bubble_sort(arr):
    """Bubble sort with operation recording"""
    operations = []
    a = arr.copy()
    n = len(a)

    for i in range(n):
        for j in range(0, n - i - 1):
            operations.append((a.copy(), [j, j + 1], 'compare', 0))
            if a[j] > a[j + 1]:
                diff = abs(a[j] - a[j + 1])
                a[j], a[j + 1] = a[j + 1], a[j]
                operations.append((a.copy(), [j, j + 1], 'swap', diff))
        operations.append((a.copy(), [n - i - 1], 'sorted', 0))

    operations.append((a.copy(), list(range(n)), 'done', 0))
    return operations

def selection_sort(arr):
    """Selection sort with operation recording"""
    operations = []
    a = arr.copy()
    n = len(a)

    for i in range(n):
        min_idx = i
        for j in range(i + 1, n):
            operations.append((a.copy(), [min_idx, j], 'compare', 0))
            if a[j] < a[min_idx]:
                min_idx = j
        if min_idx != i:
            diff = abs(a[i] - a[min_idx])
            a[i], a[min_idx] = a[min_idx], a[i]
            operations.append((a.copy(), [i, min_idx], 'swap', diff))
        operations.append((a.copy(), [i], 'sorted', 0))

    operations.append((a.copy(), list(range(n)), 'done', 0))
    return operations

def insertion_sort(arr):
    """Insertion sort with operation recording"""
    operations = []
    a = arr.copy()
    n = len(a)

    for i in range(1, n):
        key = a[i]
        j = i - 1
        while j >= 0:
            operations.append((a.copy(), [j, j + 1], 'compare', 0))
            if a[j] > key:
                diff = abs(a[j] - a[j + 1])
                a[j + 1] = a[j]
                operations.append((a.copy(), [j, j + 1], 'swap', diff))
                j -= 1
            else:
                break
        a[j + 1] = key

    operations.append((a.copy(), list(range(n)), 'done', 0))
    return operations

def quick_sort(arr):
    """Quick sort with operation recording"""
    operations = []
    a = arr.copy()

    def partition(low, high):
        pivot = a[high]
        i = low - 1

        for j in range(low, high):
            operations.append((a.copy(), [j, high], 'compare', 0))
            if a[j] <= pivot:
                i += 1
                if i != j:
                    diff = abs(a[i] - a[j])
                    a[i], a[j] = a[j], a[i]
                    operations.append((a.copy(), [i, j], 'swap', diff))

        diff = abs(a[i + 1] - a[high])
        a[i + 1], a[high] = a[high], a[i + 1]
        operations.append((a.copy(), [i + 1, high], 'swap', diff))
        return i + 1

    def quicksort_recursive(low, high):
        if low < high:
            pi = partition(low, high)
            operations.append((a.copy(), [pi], 'sorted', 0))
            quicksort_recursive(low, pi - 1)
            quicksort_recursive(pi + 1, high)

    quicksort_recursive(0, len(a) - 1)
    operations.append((a.copy(), list(range(len(a))), 'done', 0))
    return operations

def merge_sort(arr):
    """Merge sort with operation recording"""
    operations = []
    a = arr.copy()

    def merge(left, mid, right):
        left_copy = a[left:mid + 1].copy()
        right_copy = a[mid + 1:right + 1].copy()

        i = j = 0
        k = left

        while i < len(left_copy) and j < len(right_copy):
            operations.append((a.copy(), [left + i, mid + 1 + j], 'compare', 0))
            if left_copy[i] <= right_copy[j]:
                diff = abs(a[k] - left_copy[i]) if a[k] != left_copy[i] else 0
                a[k] = left_copy[i]
                i += 1
            else:
                diff = abs(a[k] - right_copy[j]) if a[k] != right_copy[j] else 0
                a[k] = right_copy[j]
                j += 1
            operations.append((a.copy(), [k], 'swap', diff))
            k += 1

        while i < len(left_copy):
            diff = abs(a[k] - left_copy[i]) if a[k] != left_copy[i] else 0
            a[k] = left_copy[i]
            operations.append((a.copy(), [k], 'swap', diff))
            i += 1
            k += 1

        while j < len(right_copy):
            diff = abs(a[k] - right_copy[j]) if a[k] != right_copy[j] else 0
            a[k] = right_copy[j]
            operations.append((a.copy(), [k], 'swap', diff))
            j += 1
            k += 1

    def mergesort_recursive(left, right):
        if left < right:
            mid = (left + right) // 2
            mergesort_recursive(left, mid)
            mergesort_recursive(mid + 1, right)
            merge(left, mid, right)

    mergesort_recursive(0, len(a) - 1)
    operations.append((a.copy(), list(range(len(a))), 'done', 0))
    return operations

# Define sorting algorithms
algorithms = [
    ("Bubble Sort", bubble_sort),
    ("Selection Sort", selection_sort),
    ("Insertion Sort", insertion_sort),
    ("Quick Sort", quick_sort),
    ("Merge Sort", merge_sort),
]

# Generate operations for all algorithms
all_operations = []
for name, func in algorithms:
    ops = func(original_data)
    all_operations.append((name, ops))
    print(f"{name}: {len(ops)} operations")

# Calculate total simulated time for each algorithm
print("\n--- Timing Summary ---")
print(f"Fixed delay per operation: {DELAY_MS} ms")
for name, ops in all_operations:
    total_time = len(ops) * DELAY_MS / 1000
    swaps = sum(1 for op in ops if op[2] == 'swap')
    print(f"{name}: {len(ops)} operations ({swaps} swaps) = {total_time:.2f} seconds")

print("\n🔊 Sound enabled: beep frequency proportional to swap difference")
print(f"   Frequency range: {MIN_FREQ} Hz (small diff) to {MAX_FREQ} Hz (large diff)")

# Now visualize each algorithm one at a time
current_algo_idx = 0
current_op_idx = 0
pause_frames = 30

fig, ax = plt.subplots(figsize=(12, 6))
bars = ax.bar(range(NUM_BARS), original_data, color='steelblue', edgecolor='black')
title = ax.set_title(f"Sorting Visualization with Sound\nStarting...", fontsize=14)
ax.set_xlabel("Index")
ax.set_ylabel("Value")
ax.set_ylim(0, 110)

# Add operation counter text
op_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=10,
                  verticalalignment='top', fontfamily='monospace')

# Color scheme
COLOR_DEFAULT = 'steelblue'
COLOR_COMPARE = 'orange'
COLOR_SWAP = 'red'
COLOR_SORTED = 'green'
COLOR_DONE = 'limegreen'

def update(frame):
    global current_algo_idx, current_op_idx, pause_frames

    if current_algo_idx >= len(all_operations):
        title.set_text("All Sorting Algorithms Complete!")
        return bars,

    name, ops = all_operations[current_algo_idx]

    if current_op_idx >= len(ops):
        if pause_frames > 0:
            pause_frames -= 1
            return bars,

        current_algo_idx += 1
        current_op_idx = 0
        pause_frames = 30

        if current_algo_idx < len(all_operations):
            for i, bar in enumerate(bars):
                bar.set_height(original_data[i])
                bar.set_color(COLOR_DEFAULT)
        return bars,

    # Get current operation
    arr_state, highlighted, op_type, swap_diff = ops[current_op_idx]

    # Play sound for swap operations
    if op_type == 'swap' and swap_diff > 0:
        freq = get_frequency_for_diff(swap_diff)
        play_beep(freq)

    # Update bar heights
    for i, bar in enumerate(bars):
        bar.set_height(arr_state[i])

        if op_type == 'done':
            bar.set_color(COLOR_DONE)
        elif i in highlighted:
            if op_type == 'compare':
                bar.set_color(COLOR_COMPARE)
            elif op_type == 'swap':
                bar.set_color(COLOR_SWAP)
            elif op_type == 'sorted':
                bar.set_color(COLOR_SORTED)
        else:
            bar.set_color(COLOR_DEFAULT)

    # Update title and operation counter
    total_ops = len(ops)
    title.set_text(f"{name}\nOperation {current_op_idx + 1} / {total_ops}")

    # Show legend
    op_text.set_text(f"Orange: Compare\nRed: Swap (with sound)\nGreen: Sorted\n\nTotal operations: {total_ops}\nTime: {total_ops * DELAY_MS / 1000:.2f}s")

    current_op_idx += 1
    return bars,

# Create animation
total_frames = sum(len(ops) + 30 for _, ops in all_operations) + 100
ani = FuncAnimation(fig, update, frames=total_frames,
                    interval=DELAY_MS, blit=False, repeat=False)

plt.tight_layout()
plt.show()
