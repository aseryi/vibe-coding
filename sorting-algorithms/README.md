# Sorting Algorithms Visualizer

A visual and auditory demonstration of different sorting algorithms.

## Features

- Visualizes 5 sorting algorithms: Bubble Sort, Selection Sort, Insertion Sort, Quick Sort, Merge Sort
- Color-coded operations: orange (compare), red (swap), green (sorted)
- Sound feedback: beep frequency proportional to the difference between swapped values
- Tracks operation count and timing for each algorithm

## Files

- `sorting_visualizer.py` - Interactive visualization with sound
- `sorting_visualizer_record.py` - Script to record the animation as MP4
- `sorting_algorithms.mp4` - Pre-recorded video of all algorithms

## Run

```bash
# Interactive with sound
python sorting_visualizer.py

# Record to video
python sorting_visualizer_record.py
```

## Requirements

- numpy
- matplotlib
- sounddevice (for audio in interactive mode)
- ffmpeg (for video recording)

## Algorithm Comparison

| Algorithm      | Operations | Swaps | Time Complexity |
|---------------|------------|-------|-----------------|
| Bubble Sort   | 324        | 113   | O(n²)           |
| Selection Sort| 226        | 15    | O(n²)           |
| Insertion Sort| 243        | 113   | O(n²)           |
| Quick Sort    | 128        | 25    | O(n log n)      |
| Merge Sort    | 155        | 88    | O(n log n)      |

*Results for 20 random elements*
