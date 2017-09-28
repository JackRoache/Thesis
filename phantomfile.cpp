#include "phantomfile.h"
#include <QFile>
#include <QDataStream>


PhantomFile::PhantomFile(QString filename, int width, int height, int frames):
    filename(filename),
    width(width),
    height(height),
    frames(frames)
{
}

QByteArray PhantomFile::getSlice(int slice)
{
    Q_ASSERT(slice < frames);
    int offset = slice * width * height;
    int frameSize = width * height;

    QFile file(filename);
    Q_ASSERT(file.open(QFile::ReadOnly));
    Q_ASSERT(file.seek(offset));


    return file.read(frameSize);
}
