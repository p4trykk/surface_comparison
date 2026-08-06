#!/usr/bin/env python3
#using in terminal: python roi_drawing.py input\file\path
import argparse
from PIL import Image, ImageDraw
import os


def draw_circle(input_path,output_path=None,center=(90, 45),radius=25,color=(0, 0, 0),width=1.5):
    img = Image.open(input_path).convert("RGB")
    draw = ImageDraw.Draw(img)
    x, y = map(int, center)
    radius = int(radius)
    width = int(width)
    bbox = [x - radius,y - radius, x + radius, y + radius]

    draw.ellipse(bbox, outline=color, width=width)
    if output_path is None:
        filename = os.path.basename(input_path)
        name, ext = os.path.splitext(filename)
        output_path = f"{name}_{ext}"


    img.save(output_path)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("input")

    parser.add_argument("-o", "--output", default=None) #using -o select output path

    parser.add_argument("-r", "--radius", type=int, default=25) #-r change circle radius 

    args = parser.parse_args()

    draw_circle(input_path=args.input,output_path=args.output, center=(90, 45),radius=args.radius)
