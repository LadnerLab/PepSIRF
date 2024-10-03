---
layout: default
title: Docker
permalink: /docker/
---

### Building the container

```
docker-compose up --build
```

**Excpected Error**:

```Error: Invalid module name entered```

This just means that the pepsirf command was not given arguments. You can ignore this since you will be providing arguments when running the image.

### Running the PepSIRF image

```
docker run --mount type=bind,src=<path/to/local/directory>,target=/app/<new_directory> \
pepsirf [ --help | module_name <module_args*> ]
```

**Note:** 
- Replace <path/to/local/directory> with the actual path to your local directory.
- Replace <new_directory> with the desired name for the directory in the container.
- If you are specifying a file (for reading or writing) from the source directory, use /app/<new_directory/file>
- Make sure that Docker is granted permission for file sharing with the local directory.
