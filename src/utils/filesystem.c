#include <stdio.h>
#include <errno.h>
#include <string.h>

#if defined(_WIN32) || defined(_WIN64)

/* -------- Windows -------- */
#include <direct.h>
#include <sys/stat.h>

int ensure_directory_exists(const char *path)
{
    struct _stat st;

    if (_stat(path, &st) == 0) {
        if (st.st_mode & _S_IFDIR)
            return 0;
        else {
            fprintf(stderr, "Path exists but is not a directory: %s\n", path);
            return -1;
        }
    }

    if (_mkdir(path) != 0) {
        perror("mkdir failed");
        return -1;
    }

    return 0;
}

#else

/* -------- POSIX (Linux, macOS) -------- */
#include <sys/stat.h>
#include <sys/types.h>

int ensure_directory_exists(const char *path)
{
    struct stat st;

    if (stat(path, &st) == 0) {
        if (S_ISDIR(st.st_mode))
            return 0;
        else {
            fprintf(stderr, "Path exists but is not a directory: %s\n", path);
            return -1;
        }
    }

    if (mkdir(path, 0755) != 0) {
        perror("mkdir failed");
        return -1;
    }

    return 0;
}

#endif
