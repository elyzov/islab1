#ifndef FEMNAMESPACE_H
#define FEMNAMESPACE_H

#include <map>
#include <string>

namespace FEM
{
/*! Типы файлов. */
    enum FileType {
        FT_SIZE,            //!< Текстовый файл.
        FT_NVTR,            //!< Двоичный файл.
        FT_NVTRSTUFF,       //!< Двоичный файл.
        FT_XY,              //!< Двоичный файл.
        FT_NVTRFST,         //!< Двоичный файл.
        FT_MU,              //!< Текстовый файл.
        FT_J,               //!< Текстовый файл.
        FT_WORKDIR,         //!< Рабочая дирректория.
        FT_RES,             //!< Двоичный файл.
        FT_B_RES            //!< Двоичный файл.
    };
    enum ActionType {
        AT_LOCAL,           //!< Добавление локальный матриц.
        AT_FIRST            //!< Учет первых краевых условий.
    };

    typedef std::map<FileType, std::string> FileMap;
}

#endif // FEMNAMESPACE_H
