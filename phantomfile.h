#ifndef PHANTOMFILE_H
#define PHANTOMFILE_H
#include <QString>

class PhantomFile
{
public:
    PhantomFile(QString filename, int width, int height, int frames);

    QByteArray getSlice(int slice);

    int getWidth(){
        return width;
    }
    int getHeight() {
        return height;
    }
    int getFrameSize() {
        return height * width;
    }

private:
    QString filename;
    int width;
    int height;
    int frames;
};

#endif // PHANTOMFILE_H
