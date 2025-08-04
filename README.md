# facade_recognition

This repository contains a Matlab library for recognising from an aerial LiDAR point cloud via a planes fitting algorithm based on the Hough transform.

## Content of the repository

The ```code``` directory contains our method that is described in Section 4.2 of the paper "PBF-FR: Partitioning Beyond Footprints for Facade Recognition in Urban Point Clouds".

The ```LAS``` directory contains the pont clouds that ... . 


## How to use it
To use this method, you can simply run the ```main.m``` file in Matlab. 

The input is a set of .las files containing the point clouds representing one or more buildings (a chunk), possibly adjacent to each other, of a city (or a part of it).
Each point cloud is processed and the method identifies the points belonging to facades.  

For each chuck, the output is a .txt file that cointains the label associated to each point of the original point cloud: 1 if the point belongs to a facade, 0 otherwise.

