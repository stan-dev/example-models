# GitHub Pages Deployment System

This directory contains tools for the automated deployment of example-models to GitHub Pages with interactive Stan Playground embeds.

## Overview

The `generate_index_files_for_ghp_deploy.py` script automatically generates `index.md` files during the GitHub Actions workflow. It has two main functions:

1. For directories containing `.stan` files:
   - Creates index files with Stan Playground embed code
   - Includes any corresponding `.data.json` file if present (optional)
   - Enables interactive model viewing/editing directly in the browser
   - Stan files without `.data.json` still get embedded with empty data

2. For other directories:
   - Creates simple index files with links to subdirectories and files
   - Maintains navigation structure through the repository

## How it Works

The script is executed automatically by GitHub Actions when changes are pushed to the repository:

1. Recursively processes all directories
2. For each directory:
   - Finds all `.stan` files and optional `.data.json` pairs
   - Generates appropriate index.md content:
     - Stan files get interactive playground embeds (with or without data)
     - Other directories get file/folder navigation links
   - Creates or updates the index.md file

This automated process enables the interactive browser experience at https://magland.github.io/example-models
