
#include <AMReX_PPMUtil.H>
#include <cstdio>
#include <cstdlib>

namespace amrex {

namespace {
    constexpr int NCOLOR = 256;
    constexpr char PGM_MAGIC1 = 'P';
//    constexpr char RPGM_MAGIC2 = '5';
    constexpr char RPGM_MAGIC3 = '6';
}

int loadPalette (const std::string& filename,
                 Array<unsigned char,NCOLOR>& r, Array<unsigned char,NCOLOR>& g,
                 Array<unsigned char,NCOLOR>& b, Array<unsigned char,NCOLOR>& a)
{
    FILE* fp = std::fopen(filename.c_str(), "rb");
    if (!fp) {
        amrex::Abort("loadPalette: cannot open "+filename);
    }

    std::fseek(fp, 0, SEEK_END);
    long length = std::ftell(fp);
    std::fseek(fp, 0, SEEK_SET);

    /* check for RGB or RGBA palette */
    int num_elements = length/(NCOLOR*sizeof(unsigned char));

    if ( num_elements != 3 && num_elements != 4 )
    {
        amrex::Abort("loadPalette: cannot process palettel file " + filename
                     + ", num(r,g,b,a) = " + std::to_string(num_elements)
                     + ", must be 3 or 4");
    }

    if (std::fread(r.data(), 1, NCOLOR, fp) != NCOLOR)
    {
        amrex::Abort("loadPalette: fread() failed to read R");
    }
    if (std::fread(g.data(), 1, NCOLOR, fp) != NCOLOR)
    {
        amrex::Abort("loadPalette: fread() failed to read G");
    }
    if (std::fread(b.data(), 1, NCOLOR, fp) != NCOLOR)
    {
        amrex::Abort("loadPalette: fread() failed to read B");
    }

    if ( num_elements == 4 ) 
    {
        if (std::fread(a.data(), 1, NCOLOR, fp) != NCOLOR)
        {
            amrex::Abort("loadPalette: fread() failed to read A");
        }
    }

    std::fclose(fp);

    return num_elements;
}

void storePPM (const std::string& filename,
               unsigned char const* data, int width, int height,
               Array<unsigned char,256> const& r,
               Array<unsigned char,256> const& g,
               Array<unsigned char,256> const& b)
{
    FILE* fp = std::fopen(filename.c_str(), "w");
    if (!fp) {
        amrex::Abort("storePPM: cannot open output file "+filename);
    }

    unsigned char* image = (unsigned char*) std::malloc(3*width*height*sizeof(unsigned char));
    if (image == nullptr) {
        amrex::Abort("storePPM: malloc failed");
    }

    for (int i = 0; i < width*height; ++i) {
        int j = static_cast<int>(data[i]);
        AMREX_ASSERT(j >= 0 && j < 256);
        image[3*i+0] = r[j];
        image[3*i+1] = g[j];
        image[3*i+2] = b[j];
    }

    std::fprintf(fp, "%c%c\n%d %d\n%d\n", PGM_MAGIC1, RPGM_MAGIC3, width, height, 255);
    std::fwrite(image, 1, 3*width*height, fp);

    std::free(image);
    std::fclose(fp);
}

}
