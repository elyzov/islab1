#include "femdatgen.h"
#include <algorithm>
#include <math.h>

FemDatGen::FemDatGen()
{
    fillFileMap();
}

FemDatGen::FemDatGen(const std::string &file_in, real zoom)
{
    fillFileMap();
    generate(file_in, zoom);
}

FemDatGen::FemDatGen(const std::string &file_in, const std::string & path_out, real zoom)
{
    fillFileMap();
    generate(file_in, path_out, zoom);
}

void FemDatGen::generate(const std::string &file_in, const std::string & path_out, real zoom)
{
    input(file_in, zoom);
    startGen();
    dataOutput(path_out);
}

void FemDatGen::generate(const std::string &file_in, real zoom)
{
    input(file_in, zoom);
    startGen();
}

void FemDatGen::dataOutput(const std::string &path_out)
{
    set_file(FT_WORKDIR, path_out);
    output();
}

void FemDatGen::clear()
{
    areas.clear();
    stuff.clear();
    gridX.clear();
    gridY.clear();
    nvtr.clear();
    xy.clear();
    nvtrFst.clear();
}

void FemDatGen::fillFileMap()
{
    files[FT_SIZE] = "inf2tr.dat";
    files[FT_NVTR] = "nvtr.dat";
    files[FT_NVTRSTUFF] = "nvkat2d.dat";
    files[FT_XY] = "rz.dat";
    files[FT_NVTRFST] = "l1.dat";
    files[FT_MU] = "mu";
    files[FT_J] = "toku";
}

void FemDatGen::input(const std::string &fname, real koef)
{
    std::ifstream fin(fname.c_str());
    if (!fin) throw Fem_OpenFile_Error(fname);
    fin.exceptions(std::ifstream::failbit |
                   std::ifstream::badbit  |
                   std::ifstream::eofbit);

    FemArea ar;
    FemGrid gr;

    fin >> kolMat;
    for (int i = 0; i < kolMat; i++)
    {
        int nmat;
        real muo;
        fin >> nmat >> muo;
        stuff[nmat] = FemStuff(muo);
    }

    fin >> nSources;
    for (int i = 0; i < nSources; i++)
    {
        real x, y, pow;
        fin >> x >> y >> pow;
        pointSrc.push_back(std::make_pair(FemPoint(x, y), pow));
    }


    fin >> x0 >> kolX;
    for (int i = 0; i < kolX; i++)
    {
        fin >> gr.coord >> gr.hMin >> gr.dh >> gr.sh;
        gridX.push_back(gr);
    }

    fin >> y0 >> kolY;
    for (int i = 0; i < kolY; i++)
    {
        fin >> gr.coord >> gr.hMin >> gr.dh >> gr.sh;
        gridY.push_back(gr);
    }

    if (koef)
    {
        for (iFemGridList it = gridX.begin(); it != gridX.end(); it++) *it *= koef;
        for (iFemGridList it = gridY.begin(); it != gridY.end(); it++) *it *= koef;
        x0 *= koef;
        y0 *= koef;
    }

    fin >> kolObl;
    for (int i = 0; i < kolObl; i++)
    {
        int nx0, nx1, ny0, ny1;
        fin >> nx0 >> nx1 >> ny0 >> ny1 >> ar.nmat;
        ar.x0 = nx0 ? gridX[nx0-1].coord : x0;
        ar.x1 = gridX[nx1-1].coord;
        ar.y0 = ny0 ? gridY[ny0-1].coord : y0;
        ar.y1 = gridY[ny1-1].coord;
        areas.push_back(ar);
    }

    fin >> doubleToX >> doubleToY;
    fin.close();
}

void FemDatGen::toTelmaGenerate(const std::string &fileIn, const std::string &fileOut, real zoom)
{
    input(fileIn, zoom);
    std::ofstream fout(fileOut.c_str());
    if (!fout) throw Fem_OpenFile_Error(fileOut);

    fout << areas.size() << '\n';
    for (size_t i = 0; i < areas.size(); i++)
    {
        fout << areas[i].x0 << '\t' << areas[i].x1 << '\t'
             << areas[i].y0 << '\t' << areas[i].y1 << '\t'
             << stuff[areas[i].nmat].sig() << '\t'
//             << stuff[areas[i].nmat].j() << '\t'
             << areas[i].nmat << '\n';
    }

/*! Out X. */
    fout << x0 << '\t' << gridX.size() << '\n';
    for (size_t i = 0; i < gridX.size(); i++)
    {
        fout << gridX[i].coord << ' ';
    }
    fout << '\n';
    for (size_t i = 0; i < gridX.size(); i++)
    {
        fout << gridX[i].hMin << ' ';
    }
    fout << '\n';
    for (size_t i = 0; i < gridX.size(); i++)
    {
        fout << gridX[i].dh << ' ';
    }
    fout << '\n';
    for (size_t i = 0; i < gridX.size(); i++)
    {
        fout << gridX[i].sh << ' ';
    }
    fout << '\n';

/*! Out Y. */
    fout << y0 * 1e-2 << '\t' << gridY.size() << '\n';
    for (size_t i = 0; i < gridY.size(); i++)
    {
        fout << gridY[i].coord << ' ';
    }
    fout << '\n';
    for (size_t i = 0; i < gridY.size(); i++)
    {
        fout << gridY[i].hMin << ' ';
    }
    fout << '\n';
    for (size_t i = 0; i < gridY.size(); i++)
    {
        fout << gridY[i].dh << ' ';
    }
    fout << '\n';
    for (size_t i = 0; i < gridY.size(); i++)
    {
        fout << gridY[i].sh << ' ';
    }
    fout << '\n';

    fout << doubleToX << ' ' << doubleToY << '\n';
}

void FemDatGen::output()
{
    std::ofstream fout;
    std::ios_base::openmode binary_mode = std::ios_base::out |
                                          std::ios_base::trunc |
                                          std::ios_base::binary;

/*! Out sizes. */
    fout.open(file(FT_SIZE).c_str());
    if (!fout) throw Fem_OpenFile_Error(file(FT_SIZE));
    fout << " islau=       0 indku1=       0 indfro=       1\n"
            " kuzlov=    " << xy.size() <<
            " ktr=    " << nvtr.size() <<
            " kt1=     " << nvtrFst.size() <<
            " kreb2=       0  kreb3=       0\n"
            " kisrr1=       2 kisrr2=       2 kisrr3=       2  kbrsr=       8\n"
             "kreb4=       0\n";
    //fout << xy.size() << " " << nvtr.size() << " " << nvtrFst.size();
    fout.close();

/*! Out nvtr. */
    fout.open(file(FT_NVTR).c_str(), binary_mode);
    if (!fout) throw Fem_OpenFile_Error(files[FT_NVTR]);
    for (int i = 0, end = nvtr.size(); i != end; i++)
    {
        long data[4];
        data[0] = (long)(nvtr[i][2] + 1);
        data[1] = (long)(nvtr[i][3] + 1);
        data[2] = (long)(nvtr[i][0] + 1);
        data[3] = (long)(nvtr[i][1] + 1);
        fout.write((char *)data, 4*sizeof(long));
        data[0] = 0; data[1] = 1;
        fout.write((char *)data, 2*sizeof(long));
    }
    fout.close();

/*! Out stuff for finite elements. */
    fout.open(file(FT_NVTRSTUFF).c_str(), binary_mode);
    if (!fout) throw Fem_OpenFile_Error(files[FT_NVTRSTUFF]);
    for (int i = 0, end = nvtr.size(); i != end; i++)
    {
        long data = nvtr[i][STUFF_IND];
        fout.write((char *)&data, sizeof(long));
    }
    fout.close();

/*! Out coords of points. */
    fout.open(file(FT_XY).c_str(), binary_mode);
    if (!fout) throw Fem_OpenFile_Error(files[FT_XY]);
    for (int i = 0, end = xy.size(); i != end; i++)
    {
        double data[2];
        data[0] = xy[i].x; data[1] = xy[i].y;
        fout.write((char *)data, 2*sizeof(double));
    }
    fout.close();

/*! Out first boundary conditions. */
    fout.open(file(FT_NVTRFST).c_str(), binary_mode);
    if (!fout) throw Fem_OpenFile_Error(files[FT_NVTRFST]);
    for (int i = 0, end = nvtrFst.size(); i != end; i++)
    {
        long data = nvtrFst[i] + 1;
        fout.write((char *)&data, sizeof(long));
    }
    fout.close();

/*! Out mu0 from stuff. */
    fout.open(file(FT_MU).c_str());
    if (!fout) throw Fem_OpenFile_Error(files[FT_MU]);
    for (iFemStuffList it = stuff.begin(); it != stuff.end(); ++it)
    {
        fout << it->first << ' ' << it->second.sig() << '\n';
    }
    fout.close();

///*! Out J from stuff. */
//    fout.open(file(FT_J).c_str());
//    if (!fout) throw Fem_OpenFile_Error(files[FT_J]);
//    for (iFemStuffList it = stuff.begin(); it != stuff.end(); ++it)
//    {
//        fout << it->first << ' ' << it->second.j() << '\n';
//    }
//    fout.close();
}

void FemDatGen::startGen()
{
    doubleTo(gridX, doubleToX);
    doubleTo(gridY, doubleToY);
    divAxis(x0, gridX, vx);
    divAxis(y0, gridY, vy);

    ndx = vx.size();    //!< number of nodes on OX axis.
    nx  = ndx - 1;      //!< number of segments on OX axis.
    ndy = vy.size();    //!< number of nodes on OY axis.
    ny  = ndy - 1;      //!< number of segments on OY axis.

//    std::sort(areas.begin(), areas.end());

    intVector elem(NVTR_DIM);
    int nArea = 0;          //!< number of current area.
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            FemPoint p(vx[i], vy[j]);
            if (!pointInArea(p, nArea)) nArea = findArea(p);

            int gl = i + ndx * j;
            elem[0] = gl;
            elem[1] = gl + 1;
            elem[2] = gl + ndx;
            elem[3] = gl + ndx + 1;
            elem[4] = areas[nArea].nmat;

            xy.push_back(p);
            nvtr.push_back(elem);
        }
        xy.push_back(FemPoint(vx[nx], vy[j]));
    }
    for (int i = 0; i < ndx; i++)
    {
        xy.push_back(FemPoint(vx[i], vy[ny]));
    }

/*! Calculate index of point source. */
    for (int i = 0; i < nSources; i++)
    {
        sources.push_back(std::make_pair(
                              getPointInd(pointSrc[i].first),
                              pointSrc[i].second)
                          );
    }
     //!< Всегда находится в левом верхнем углу.

/*! Generating array of first boundary conditions. */
    /*! Нижняя граница. */
    for (int i = 0; i < ndx; i++) nvtrFst.push_back(i);
    for (int j = 0; j < ndy; j++)
    {
//        nvtrFst.push_back(j * ndx); //!< Левая граница.
        nvtrFst.push_back(j * ndx + nx); //!< Правая граница.
    }
    /*! Верхняя граница. */
//    for (int i = ny * ndx, end = ndy * ndx; i < end; i++)
//        nvtrFst.push_back(i);
}

int FemDatGen::divAxis(real startCoord, const FemGridList &grid, rvector & v)
{
    v.clear();
    v.push_back(startCoord);
    for (ciFemGridList it = grid.begin(); it != grid.end(); it++)
    {
        real h = it->hMin;
        if (it->dh < 1) throw Fem_Grid_Error("dh < 1");
        if (it->sh < 0)
        {
            real coord = it->coord;
            rlist lv;
            while (startCoord + cmp_eps < coord)
            {
                lv.push_front(coord);
                coord += h;
                h *= it->dh;
            }
            v.insert(v.end(), lv.begin(), lv.end());
        }
        else
        {
            real coord = startCoord + h;
            while (coord + cmp_eps < it->coord)
            {
                v.push_back(coord);
                coord += h;
                h *= it->dh;
            }
            v.push_back(it->coord);
            startCoord = it->coord;
        }
    }
    return 0;
}

bool FemDatGen::pointInArea(const FemPoint &p, int nArea)
{
    const FemArea & ar = areas[nArea];
    if (p.x >= ar.x0 && p.x < ar.x1
     && p.y >= ar.y0 && p.y < ar.y1) return true;
    else return false;
}

void FemDatGen::doubleTo(FemGridList & grid, int times)
{
    if (times <= 0) return;
    for (iFemGridList it = grid.begin(); it != grid.end(); it++)
    {
        if (fabs(it->dh - 1.) < cmp_eps)
        {
            it->hMin /= times * 2.;
        }
        else
        {
            it->dh = pow(it->dh, 1. / times);
            it->hMin /= it->dh;
        }
    }
}

size_t FemDatGen::findArea(const FemPoint &p)
{
    for (size_t i = 0, end = areas.size(); i < end; i++)
    {
        if (pointInArea(p, i)) return i;
    }
    throw Fem_Grid_Error("Area not found.");
}

void FemDatGen::set_file(FileType ftype, const std::string & fname)
{
    files[ftype] = fname;
    if (ftype == FT_WORKDIR) files[ftype] = correct_path(fname);
}

std::string FemDatGen::correct_path(const std::string & path)
{
    std::string res = path;
    size_t end = res.size() - 1;
    if (res[end] != '/') res.push_back('/');
    return res;
}

inline const std::string FemDatGen::file(FileType ftype)
{
    return files[FT_WORKDIR] + files[ftype];
}

int FemDatGen::getPointInd(const FemPoint &p)
{
    int npx, npy;
    irvector ix = std::find(vx.begin(), vx.end(), p.x);
    if (ix == vx.end()) throw Fem_IO_Error("Point source must be located at grid");
    else npx = (int)(ix - vx.begin());
    irvector iy = std::find(vy.begin(), vy.end(), p.y);
    if (iy == vy.end()) throw Fem_IO_Error("Point source must be located at grid");
    else npy = (int)(iy - vy.begin());
    return npy * ndx + npx;
}
