# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

# Configuration
NUM_BARS = 20
DELAY_MS = 80  # Delay between operations
FPS = 1000 // DELAY_MS  # Convert to frames per second

# Generate random data
np.random.seed(42)
original_data = np.random.randint(1, 100, NUM_BARS)
max_diff = 100

# Store all operations for each sorting algorithm
def bubble_sort(arr):
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

print("\n--- Recording Video ---")
print(f"FPS: {FPS}")

# Flatten all operations into a single sequence with algorithm info
frames_data = []
for name, ops in all_operations:
    # Add pause frames at start showing algorithm name
    for _ in range(15):
        frames_data.append((name, original_data.copy(), [], 'start', 0, 0, len(ops)))

    for idx, (arr_state, highlighted, op_type, swap_diff) in enumerate(ops):
        frames_data.append((name, arr_state, highlighted, op_type, swap_diff, idx + 1, len(ops)))

    # Add pause frames at end
    for _ in range(15):
        frames_data.append((name, ops[-1][0], list(range(NUM_BARS)), 'done', 0, len(ops), len(ops)))

# Add final frames
for _ in range(30):
    frames_data.append(("Complete!", ops[-1][0], list(range(NUM_BARS)), 'final', 0, 0, 0))

total_frames = len(frames_data)
print(f"Total frames: {total_frames}")
print(f"Video duration: {total_frames / FPS:.1f} seconds")

# Set up the figure
fig, ax = plt.subplots(figsize=(12, 6))
bars = ax.bar(range(NUM_BARS), original_data, color='steelblue', edgecolor='black')
title = ax.set_title("Sorting Algorithms Visualization", fontsize=14)
ax.set_xlabel("Index")
ax.set_ylabel("Value")
ax.set_ylim(0, 110)

op_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=10,
                  verticalalignment='top', fontfamily='monospace')

# Color scheme
COLOR_DEFAULT = 'steelblue'
COLOR_COMPARE = 'orange'
COLOR_SWAP = 'red'
COLOR_SORTED = 'green'
COLOR_DONE = 'limegreen'

def update(frame):
    name, arr_state, highlighted, op_type, swap_diff, current_op, total_ops = frames_data[frame]

    # Update bar heights and colors
    for i, bar in enumerate(bars):
        bar.set_height(arr_state[i])

        if op_type in ['done', 'final']:
            bar.set_color(COLOR_DONE)
        elif op_type == 'start':
            bar.set_color(COLOR_DEFAULT)
        elif i in highlighted:
            if op_type == 'compare':
                bar.set_color(COLOR_COMPARE)
            elif op_type == 'swap':
                bar.set_color(COLOR_SWAP)
            elif op_type == 'sorted':
                bar.set_color(COLOR_SORTED)
        else:
            bar.set_color(COLOR_DEFAULT)

    # Update title
    if op_type == 'final':
        title.set_text("All Sorting Algorithms Complete!")
    elif op_type == 'start':
        title.set_text(f"{name}\nStarting...")
    else:
        title.set_text(f"{name}\nOperation {current_op} / {total_ops}")

    # Update legend
    if op_type != 'final':
        op_text.set_text(f"Orange: Compare\nRed: Swap\nGreen: Sorted\n\nOperations: {total_ops}")
    else:
        op_text.set_text("Done!")

    return bars,

print("\nRendering video...")
ani = FuncAnimation(fig, update, frames=total_frames, blit=False, repeat=False)

# Save as MP4
output_path = "/Users/seryi/Library/CloudStorage/GoogleDrive-andrei.seryi@gmail.com/My Drive/Claude/sorting_algorithms.mp4"
writer = FFMpegWriter(fps=FPS, metadata=dict(title='Sorting Algorithms'), bitrate=2000)
ani.save(output_path, writer=writer)

print(f"\nVideo saved to: {output_path}")
plt.close()
