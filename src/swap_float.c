float swap_float(float val)
{
    union {
        float f;
        unsigned char b[4];
    } src, dst;

    src.f = val;
    dst.b[0] = src.b[3];
    dst.b[1] = src.b[2];
    dst.b[2] = src.b[1];
    dst.b[3] = src.b[0];

    return dst.f;
}