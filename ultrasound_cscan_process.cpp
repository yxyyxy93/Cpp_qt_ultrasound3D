#include "ultrasound_cscan_process.h"
#include "qcustomplot.h"
#include "utils.h"

#include <QtWidgets>
#include <QHBoxLayout>
#include <QScrollBar>
#include <QDebug>
#include <QList>
#include <QTreeView>
#include <QStandardItemModel>
#include <algorithm>

#include "npy.hpp"

ultrasound_Cscan_process::ultrasound_Cscan_process(QWidget *parent,
                                                   QString fn,
                                                   int fs,
                                                   double fx,
                                                   double fy)
    : QWidget(parent)
    , fn(fn)
    , fs(fs)
    , fx(fx)
    , fy(fy)
    , customPlot_Ascan(nullptr)
{
    this->layout = new QVBoxLayout(this);

    // create a QMainWindow and set it as the central widget
    QMainWindow *mainWindow = new QMainWindow(this);
    mainWindow->setCentralWidget(new QWidget(mainWindow));
    this->layout->addWidget(mainWindow);

    // ************** create the push buttons and correponding labels
    this->myButton_load = new QPushButton(tr("Load data"), parent);
    layout->addWidget(this->myButton_load);
    connect(this->myButton_load,
            &QPushButton::clicked,
            this,
            &ultrasound_Cscan_process::handleButton_load);
    //
    this->myButton_save = new QPushButton(tr("Save"), this);
    connect(this->myButton_save,
            &QPushButton::clicked,
            this,
            &ultrasound_Cscan_process::handleButton_save);
    layout->addWidget(this->myButton_save);
    //
    this->myButton_loadraw = new QPushButton(tr("Load raw"), this);
    connect(this->myButton_loadraw,
            &QPushButton::clicked,
            this,
            &ultrasound_Cscan_process::handleButton_loadraw);
    layout->addWidget(this->myButton_loadraw);
    //
    this->myButton_AS = new QPushButton("calculate the analytic-signal", this);
    connect(this->myButton_AS,
            &QPushButton::clicked,
            this,
            &ultrasound_Cscan_process::handleButton_AS);
    layout->addWidget(this->myButton_AS);
    //
    this->myButton_orthoslice = new QPushButton(tr("Orthoslice"), this);
    connect(this->myButton_orthoslice,
            &QPushButton::clicked,
            this,
            &ultrasound_Cscan_process::handleButton_orthoslice);
    layout->addWidget(this->myButton_orthoslice);
    //
    this->myButton_surface = new QPushButton("Determine the surface", this);
    connect(this->myButton_surface,
            &QPushButton::clicked,
            this,
            &ultrasound_Cscan_process::handleButton_surface);
    layout->addWidget(this->myButton_surface);

    // ************** progress bar
    this->m_progressBar = new QProgressBar;
    this->m_progressBar->setRange(0, 100);
    layout->addWidget(this->m_progressBar);
    setLayout(layout);
}

ultrasound_Cscan_process::~ultrasound_Cscan_process(){
    // free any resources that were allocated by the class
    delete myButton_load;
    delete myButton_save;
    delete myButton_loadraw;
    delete myButton_orthoslice;
    delete myButton_surface;

    delete m_progressBar;

    delete myButton_alignsurface;
    delete myButton_AS;

    delete customPlot1;
    delete customPlot2;
    delete customPlot3;
    delete scrollBarX;
    delete scrollBarY;
    delete scrollBarZ;

    delete customPlot_Ascan;

    delete customPlot_frontI;
    delete customPlot_frontV;

    delete layout;
}

// ***********　ｒｅａｄ　ｄａｔａ　－　ｂｕｔｔｏｎ
void ultrasound_Cscan_process::handleButton_load()
{
    // get the file name
    this->fn = QFileDialog::getOpenFileName(nullptr,
                                            "Open file",
                                            "",
                                            "Text files (*.txt; *.tdms; *npy)");
    if (this->fn.isEmpty()) {
        qWarning() << "No file selected.";
    }
    // Determine the file type based on its extension
    QFileInfo fileInfo(this->fn);
    QString extension = fileInfo.suffix();
    QVector<QVector<double>> B_scan_double;
    QVector<double> A_scan_double;
    this->C_scan_double.clear();
    if (extension == "txt") {
        // The file is a text file
        // read the data
        QFile file(this->fn);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            qWarning() << "Could not open file:" << file.errorString();
        }
        QString contents = file.readAll().constData();
        file.close();
        // split the raw data
        QStringList contents_split;
        contents_split << contents.split("\n");
        // *************** read data
        QStringList lineData;
        // count the signals
        int count(1); // avoid divied by zero
        int lengthThreshold = 500;
        for (const QString &str : contents_split) {
            if (str.length() > lengthThreshold) {
                count++;
            }
        }
        qDebug() << "count of signals:" << count-1;
        //
        int count_cur(1);
        for (const QString &thisline: contents_split){
            // flag of the start of a B-scan
            if (thisline.length()>=5 && thisline.contains("GROUP")){
                C_scan_double.push_back(B_scan_double);
                B_scan_double.clear();
            }
            lineData.clear();
            lineData << thisline.split("\t");
            if (lineData.size() > lengthThreshold){
                count_cur++;
                m_progressBar->setValue(100 * count_cur / count);
                A_scan_double.clear();
                for (const QString &elem: lineData)
                {
                    A_scan_double.push_back(elem.toDouble());
                    // Update the progress bar
                    QCoreApplication::processEvents(); // Allow GUI updates
                }
                B_scan_double.push_back(A_scan_double);
            }
        }
        // last B_scan
        this->C_scan_double.removeFirst(); // remove the first one
        this->C_scan_double.push_back(B_scan_double);
        m_progressBar->setValue(100);
    } else if (extension == "tdms"){
        ;
    } else if (extension == "npy"){
        std::vector<unsigned long> shape {};
        bool fortran_order;
        std::vector<double> data;
        const std::string path = this->fn.toStdString();
        npy::LoadArrayFromNumpy(path, shape, fortran_order, data);
        // Define the dimensions of the output QVector !! to be improved !!!
        int x = 146;
        int y = 150;
        int z = data.size()/x/y;
        // Copy the input data into the output QVector
        int idx = 0;
        for (int i = 0; i < x; i++) {
            B_scan_double.clear();
            for (int j = 0; j < y; j++) {
                A_scan_double.clear();
                for (int k = 0; k < z; k++) {
                    A_scan_double.push_back(data[idx]);
                    idx++;
                }
                B_scan_double.push_back(A_scan_double);
            }
            m_progressBar->setValue(100 * i / y);
            // Update the progress bar
            QCoreApplication::processEvents(); // Allow GUI updates
            this->C_scan_double.push_back(B_scan_double);
        }
        // compensate to make it square
        for (int i = 0; i< y-x; i++)
            this->C_scan_double.push_back(B_scan_double);
        m_progressBar->setValue(100);
    } else {
        qDebug() << "The file type is unknown";
    }
    qDebug() << this->C_scan_double.size();
    qDebug() << this->C_scan_double[0].size();
    qDebug() << this->C_scan_double[0][0].size();
}

// *********** save and read data - button
void ultrasound_Cscan_process::handleButton_save(){
    QString filename = QFileDialog::getSaveFileName(this,
                                                    tr("Save Vector"),
                                                    QDir::homePath(),
                                                    tr("Text Files (*.txt)"));
    QFile file(filename);
    if (file.open(QIODevice::WriteOnly)) {
        QDataStream out(&file); // Create a QDataStream to write to the file

        // Write the dimensions as a header
        qint32 n1 = this->C_scan_double.size();
        qint32 n2 = this->C_scan_double[0].size();
        qint32 n3 = this->C_scan_double[0][0].size();
        out << n1 << n2 << n3;

        // Serialize the data to the stream
        out << this->C_scan_double;

        file.close();
    }
}

void ultrasound_Cscan_process::handleButton_loadraw() {
    //    QString filename = QFileDialog::getOpenFileName(this,
    //                                                    tr("Open Data"),
    //                                                    QDir::homePath(),
    //                                                    tr("Text Files (*.txt)"));
    //    qDebug() << filename;
    //    define the filename directly
    QString filename = "C:/Users/xiayang/OneDrive - UGent/lab/experiment/Woven_Samples/2-ceh206-8-p21-0/preprocessed.txt";
    if (!filename.isEmpty()) {
        QFile file(filename);
        QVector<QVector<QVector<double>>> myData;
        if (file.open(QIODevice::ReadOnly)) {
            QDataStream in(&file); // Create a QDataStream to read from the file
            // Read the dimensions from the header
            qint32 n1, n2, n3;
            in >> n1 >> n2 >> n3;
            myData.resize(n1);
            for (int i = 0; i < n1; ++i) {
                myData[i].resize(n2);
                for (int j = 0; j < n2; ++j) {
                    myData[i][j].resize(n3);
                }
            }
            // Deserialize the data from the stream
            in >> myData;
            file.close();
        }
        // asign the data
        this->C_scan_double = myData;
    }
    // done
    m_progressBar->setValue(100);
    qDebug() << this->C_scan_double.size();
    qDebug() << this->C_scan_double[0].size();
    qDebug() << this->C_scan_double[0][0].size();
    // calculate Analytic-signal
    QVector<std::complex<double>> Ascan_as;
    QVector<QVector<std::complex<double>>> Bscan_as;

    for(int i = 0; i < this->C_scan_double.size(); i++) {
        for(int j = 0; j < this->C_scan_double[i].size(); j++) {
            QVector<double> Ascan = this->C_scan_double[i][j];
            Ascan_as = analyticSignal(Ascan);
            Bscan_as.push_back(Ascan_as);
        }
        this->m_progressBar->setValue(100 * i / this->C_scan_double.size());
        // Update the progress bar
        QCoreApplication::processEvents(); // Allow GUI updates
        this->C_scan_AS.push_back(Bscan_as);
        Bscan_as.clear();
    }
    this->m_progressBar->setValue(100);
    this->myButton_AS ->setText("Analytic-signal ready!");
}

void ultrasound_Cscan_process::handleButton_AS(){
    QVector<std::complex<double>> Ascan_as;
    QVector<QVector<std::complex<double>>> Bscan_as;

    for(int i = 0; i < this->C_scan_double.size(); i++) {
        for(int j = 0; j < this->C_scan_double[i].size(); j++) {
            QVector<double> Ascan = this->C_scan_double[i][j];
            Ascan_as = analyticSignal(Ascan);
            Bscan_as.push_back(Ascan_as);
        }
        this->m_progressBar->setValue(100 * i / this->C_scan_double.size());
        // Update the progress bar
        QCoreApplication::processEvents(); // Allow GUI updates
        this->C_scan_AS.push_back(Bscan_as);
        Bscan_as.clear();
    }
    this->m_progressBar->setValue(100);
    this->myButton_AS ->setText("Analytic-signal ready!");
}

// **************** define surface
// Slot to perform the task
void ultrasound_Cscan_process::handleButton_surface()
{
    // loop through each element and compare to max value
    QVector<double> Front_surface_idx_i;
    QVector<double> Front_surface_val_i;
    for(int i = 0; i < this->C_scan_double.size(); i++){
        for(int j = 0; j < this->C_scan_double[i].size(); j++){
            // Find the maximum value and its index
            QVector<std::complex<double>> Ascan_as = this->C_scan_AS[i][j];
            // find the max
            QVector<std::complex<double>>::iterator maxElementIndex;
            maxElementIndex = std::max_element(Ascan_as.begin(),
                                               Ascan_as.end(),
                                               [](std::complex<double> a, std::complex<double> b) {
                    return std::abs(a) < std::abs(b);
        });
            Front_surface_idx_i.push_back(std::distance(Ascan_as.begin(), maxElementIndex));
            Front_surface_val_i.push_back(std::abs(*maxElementIndex));
            this->m_progressBar->setValue(100 * i / this->C_scan_double.size());
            // Update the progress bar
            QCoreApplication::processEvents(); // Allow GUI updates
        }
        // add to 2D QVector
        this->Front_surface_idx.push_back(Front_surface_idx_i);
        this->Front_surface_val.push_back(Front_surface_val_i);
        Front_surface_idx_i.clear();
        Front_surface_val_i.clear();
    }
    this->m_progressBar->setValue(100);
    // add a Qlabel
    QLabel *statusLabel_surface = new QLabel(this);
    this->layout->addWidget(statusLabel_surface);
    // When the task is finished, update the label
    statusLabel_surface->setText("Surface determination completed successfully");
    // add a button to plot the surface
    // ************** create the push buttons and correponding labels
    QPushButton* myButton_plotsurface = new QPushButton(tr("Plot front surface"), this);
    layout->addWidget(myButton_plotsurface);
    connect(myButton_plotsurface,
            &QPushButton::clicked, this,
            &ultrasound_Cscan_process::handleButton_plotsurface);
    this->myButton_alignsurface = new QPushButton(tr("Align front surface"), this);
    layout->addWidget(this->myButton_alignsurface);
    connect(this->myButton_alignsurface,
            &QPushButton::clicked, this,
            &ultrasound_Cscan_process::handleButton_alignsurface);
}

void ultrasound_Cscan_process::handleButton_alignsurface(){
    // get min of the front surface index
    int min_idx = this->C_scan_AS[0][0].size();
    for(int i = 0; i < this->Front_surface_idx.size(); i++) {
        for(int j = 0; j < this->Front_surface_idx[i].size(); j++) {
            min_idx = (min_idx>=this->Front_surface_idx[i][j])?
                        this->Front_surface_idx[i][j]:min_idx;
        }
    }
    min_idx = 100; // manual setting !!!!!!!!
    qDebug() << min_idx;
    // shift
    for(int i = 0; i < this->Front_surface_idx.size(); i++){
        for(int j = 0; j < this->Front_surface_idx[i].size(); j++){
            int front_idx = this->Front_surface_idx[i][j];
            shiftVector_1D(this->C_scan_AS[i][j],
                           front_idx-min_idx);
        }
    }
    this->myButton_alignsurface ->setText("Surface aligned !");
}

void ultrasound_Cscan_process::handleButton_plotsurface(){
    // plot the surface
    QCustomPlot *customPlot_fsurface_idx = new QCustomPlot();
    QCustomPlot *customPlot_fsurface_val = new QCustomPlot();
    //
    this->layout->addWidget(customPlot_fsurface_idx);
    this->layout->addWidget(customPlot_fsurface_val);
    // Create QCPColorMap objects
    QCPColorMap *map1 = new QCPColorMap(customPlot_fsurface_idx->xAxis,
                                        customPlot_fsurface_idx->yAxis);
    QCPColorMap *map2 = new QCPColorMap(customPlot_fsurface_val->xAxis,
                                        customPlot_fsurface_val->yAxis);
    // set the data for each QCPColorMap
    int x_size = this->Front_surface_idx.size();
    int y_size = this->Front_surface_idx[0].size();
    //
    map1->data()->setSize(x_size, y_size);
    map1->data()->setRange(QCPRange(0, x_size),
                           QCPRange(0, y_size));
    map2->data()->setSize(x_size, y_size);
    map2->data()->setRange(QCPRange(0, x_size),
                           QCPRange(0, y_size));
    //
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j) {
            map1->data()->setCell(i, j, this->Front_surface_idx[i][j]);
        }
    }
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j) {
            map2->data()->setCell(i, j, this->Front_surface_val[i][j]);
        }
    }
    map1->setGradient(QCPColorGradient::gpHot);
    map2->setGradient(QCPColorGradient::gpHot);
    map1->setInterpolate(true);
    map2->setInterpolate(true);
    // Add a color scale to the custom plot widget
    QCPColorScale *colorScale;
    colorScale = new QCPColorScale(customPlot_fsurface_idx);
    map1->setColorScale(colorScale);
    colorScale = new QCPColorScale(customPlot_fsurface_val);
    map2->setColorScale(colorScale);
    // Rescale the color map data range to fit the new data
    map1->rescaleDataRange();
    map2->rescaleDataRange();
    //
    customPlot_fsurface_idx->xAxis->setRange(0, x_size);
    customPlot_fsurface_idx->yAxis->setRange(0, y_size);
    customPlot_fsurface_val->xAxis->setRange(0, x_size);
    customPlot_fsurface_val->yAxis->setRange(0, y_size);
    // Call replot() to update the plot with the new data
    customPlot_fsurface_idx->replot();
    customPlot_fsurface_val->replot();
}

// *************ｖｉｓｕａｌｉｚａｔｉｏｎ
void ultrasound_Cscan_process::handleButton_orthoslice() {
    // ****************** create the orthoslice visual
    // Create the QCustomPlot widget
    this->customPlot1 = new QCustomPlot();
    this->customPlot2 = new QCustomPlot();
    this->customPlot3 = new QCustomPlot();
    // Qcustomplot
    this->scrollBarX = new QScrollBar();
    this->scrollBarY = new QScrollBar();
    this->scrollBarZ = new QScrollBar();
    // Set the positions and sizes of the QCustomPlot and QScrollBars
    this->customPlot1->setGeometry(50, 50, 600, 400);
    //
    this->scrollBarX->setOrientation(Qt::Horizontal);
    this->scrollBarY->setOrientation(Qt::Horizontal);
    this->scrollBarZ->setOrientation(Qt::Horizontal);
    this->scrollBarX->setGeometry(50, 470, 600, 20);
    this->scrollBarY->setGeometry(50, 500, 600, 20);
    this->scrollBarZ->setGeometry(50, 530, 600, 20);
    // Create a horizontal layout for the main window
    QWidget *rightPanel = new QWidget();
    QHBoxLayout *hLayout = new QHBoxLayout();
    rightPanel->setLayout(hLayout);
    // Create a vertical layout for 1
    QVBoxLayout *plot1l = new QVBoxLayout();
    QWidget *plot1w = new QWidget();
    plot1w->setLayout(plot1l);
    // Add some widgets to the right panel
    plot1l->addWidget(this->customPlot1);
    plot1l->addWidget(this->scrollBarX);
    // Create a vertical layout for 2
    QVBoxLayout *plot2l = new QVBoxLayout();
    QWidget *plot2w = new QWidget();
    plot2w->setLayout(plot2l);
    // Add some widgets to the right panel
    plot2l->addWidget(this->customPlot2);
    plot2l->addWidget(this->scrollBarY);
    // Create a vertical layout for 3
    QVBoxLayout *plot3l = new QVBoxLayout();
    QWidget *plot3w = new QWidget();
    plot3w->setLayout(plot3l);
    // Add some widgets to the right panel
    plot3l->addWidget(this->customPlot3);
    plot3l->addWidget(this->scrollBarZ);
    //
    hLayout->addWidget(plot1w);
    hLayout->addWidget(plot2w);
    hLayout->addWidget(plot3w);
    //
    this->layout->addWidget(rightPanel);
    // Connect the valueChanged() signals of the QScrollBars to update the plot data
    QObject::connect(this->scrollBarX, &QScrollBar::valueChanged, this,
                     &ultrasound_Cscan_process::updatePlot);
    QObject::connect(this->scrollBarY, &QScrollBar::valueChanged, this,
                     &ultrasound_Cscan_process::updatePlot);
    QObject::connect(this->scrollBarZ, &QScrollBar::valueChanged, this,
                     &ultrasound_Cscan_process::updatePlot);
    // Set the range of the QScrollBars based on the size of the data
    this->scrollBarX->setRange(0, this->C_scan_double.size() - 1);
    this->scrollBarY->setRange(0, this->C_scan_double[0].size() - 1);
    this->scrollBarZ->setRange(0, this->C_scan_double[0][0].size() - 1);
    // Set the initial values of the QScrollBars
    this->scrollBarX->setValue(this->C_scan_double.size() / 2);
    this->scrollBarY->setValue(this->C_scan_double[0].size() / 2);
    this->scrollBarZ->setValue(this->C_scan_double[0][0].size() / 2);
    // Update the plot
    updatePlot();
}

void ultrasound_Cscan_process::updatePlot() {
    int sliceX = this->scrollBarX->value();
    int sliceY = this->scrollBarY->value();
    int sliceZ = this->scrollBarZ->value();
    qDebug() << sliceX;
    qDebug() << sliceY;
    qDebug() << sliceZ;
    // Create QCPColorMap objects
    QCPColorMap *map1 = new QCPColorMap(this->customPlot1->xAxis,
                                        this->customPlot1->yAxis);
    QCPColorMap *map2 = new QCPColorMap(this->customPlot2->xAxis,
                                        this->customPlot2->yAxis);
    QCPColorMap *map3 = new QCPColorMap(this->customPlot3->xAxis,
                                        this->customPlot3->yAxis);
    // set the data for each QCPColorMap
    int x_size = this->C_scan_AS.size();
    int y_size = this->C_scan_AS[0].size();
    int z_size = this->C_scan_AS[0][0].size();
    //
    map1->data()->setSize(y_size, z_size);
    map1->data()->setRange(QCPRange(0, y_size),
                           QCPRange(0, z_size));
    map2->data()->setSize(x_size, z_size);
    map2->data()->setRange(QCPRange(0, x_size),
                           QCPRange(0, z_size));
    map3->data()->setSize(x_size, y_size);
    map3->data()->setRange(QCPRange(0, x_size),
                           QCPRange(0, y_size));
    //
    for (int i = 0; i < y_size; ++i) {
        for (int j = 0; j < z_size; ++j) {
            map1->data()->setCell(i, j, std::abs(this->C_scan_AS[sliceX][i][j]));
        }
    }

    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < z_size; ++j) {
            map2->data()->setCell(i, j, std::abs(this->C_scan_AS[i][sliceY][j]));
        }
    }

    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j) {
            map3->data()->setCell(i, j, std::abs(this->C_scan_AS[i][j][sliceZ]));
        }
    }
    //
    map3->setGradient(QCPColorGradient::gpJet);
    //
    map1->setInterpolate(true);
    map2->setInterpolate(true);
    map3->setInterpolate(true);
    // Add a color scale to the custom plot widget
    QCPColorScale *colorScale = new QCPColorScale(this->customPlot3);
    // Set the color map for the color scale
    colorScale->setDataRange(map3->dataRange());
    colorScale->setGradient(map3->gradient());
    map1->setColorScale(colorScale);
    map2->setColorScale(colorScale);
    map3->setColorScale(colorScale);
    // add a color scale:
    this->customPlot3->plotLayout()->addElement(0, 1, colorScale); // add it to the right of the main axis rect
    colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)// associate the color map with the color scale
    colorScale->axis()->setLabel("Amp. (arb.)");
    // Rescale the color map data range to fit the new data
    map1->rescaleDataRange(true);
    map2->rescaleDataRange(true);
    map3->rescaleDataRange(true);
    //
    map1->rescaleAxes();
    map2->rescaleAxes();
    map3->rescaleAxes();
    //
    this->customPlot1->xAxis->setRange(0, y_size);
    this->customPlot1->yAxis->setRange(0, z_size);
    this->customPlot2->xAxis->setRange(0, x_size);
    this->customPlot2->yAxis->setRange(0, z_size);
    this->customPlot3->xAxis->setRange(0, x_size);
    this->customPlot3->yAxis->setRange(0, y_size);
    // Call replot() to update the plot with the new data
    this->customPlot1->replot();
    this->customPlot2->replot();
    this->customPlot3->replot();
    // Connect the plot's mouse press signal to the onCustomPlotClicked slot
    connect(this->customPlot3, &QCustomPlot::mousePress,
            this, &ultrasound_Cscan_process::onCustomPlotClicked_Cscan);
}

void ultrasound_Cscan_process::onCustomPlotClicked_Cscan(QMouseEvent* event)
{
    static QElapsedTimer lastClickTime;
    int minClickInterval = 500; // minimum time between clicks in milliseconds
    //
    if (event->button() == Qt::LeftButton)
    {
        if (lastClickTime.isValid() && lastClickTime.elapsed() < minClickInterval)
        {
            // ignore click if it's too soon after the last one
            return;
        }
        lastClickTime.start();
        // handle the click event
        int x = this->customPlot3->xAxis->pixelToCoord(event->pos().x());
        int y = this->customPlot3->yAxis->pixelToCoord(event->pos().y());
        qDebug() << "Clicked at (" << x << "," << y << ")";
        // plot Ascan
        QVector<double> signal = this->C_scan_double[x][y];
        QVector<double> time;
        for (int i = 0; i < signal.size(); ++i) {
            time.append(i);
        }
        // calculate the analytic-signal
        QVector<std::complex<double>> Ascan_as;
        Ascan_as = analyticSignal(signal);
        QVector<double> Ascan_as_abs;
        for (const auto& element : Ascan_as) {
            double absoluteValue = std::abs(element);
            Ascan_as_abs.append(absoluteValue);
        }
        // create a new
        if (this->customPlot_Ascan==nullptr){
            this->customPlot_Ascan = new QCustomPlot();
            this->layout->addWidget(this->customPlot_Ascan);
        }
        else
            this->customPlot_Ascan->clearGraphs();
        // add two new graphs and set their look:
        this->customPlot_Ascan->addGraph();
        this->customPlot_Ascan->graph(0)->setPen(QPen(Qt::blue)); // line color blue for first graph
        this->customPlot_Ascan->graph(0)->setBrush(QBrush(QColor(0, 0, 255, 20))); // first graph will be filled with translucent blue
        this->customPlot_Ascan->addGraph();
        this->customPlot_Ascan->graph(1)->setPen(QPen(Qt::red)); // line color red for second graph
        // configure right and top axis to show ticks but no labels:
        // (see QCPAxisRect::setupFullAxesBox for a quicker method to do this)
        this->customPlot_Ascan->xAxis2->setVisible(true);
        this->customPlot_Ascan->xAxis2->setTickLabels(false);
        this->customPlot_Ascan->yAxis2->setVisible(true);
        this->customPlot_Ascan->yAxis2->setTickLabels(false);
        // make left and bottom axes always transfer their ranges to right and top axes:
        connect(this->customPlot_Ascan->xAxis,
                SIGNAL(rangeChanged(QCPRange)),
                this->customPlot_Ascan->xAxis2,
                SLOT(setRange(QCPRange)));
        connect(this->customPlot_Ascan->yAxis,
                SIGNAL(rangeChanged(QCPRange)),
                this->customPlot_Ascan->yAxis2,
                SLOT(setRange(QCPRange)));
        // pass data points to graphs:
        this->customPlot_Ascan->graph(0)->setData(time,
                                                  signal);
        this->customPlot_Ascan->graph(1)->setData(time,
                                                  Ascan_as_abs);
        // let the ranges scale themselves so graph 0 fits perfectly in the visible area:
        this->customPlot_Ascan->graph(0)->rescaleAxes();
        // same thing for graph 1, but only enlarge ranges (in case graph 1 is smaller than graph 0):
        this->customPlot_Ascan->graph(1)->rescaleAxes(true);
        // Note: we could have also just called customPlot->rescaleAxes(); instead

        // Allow user to drag axis ranges with mouse, zoom with mouse wheel and select graphs by clicking:
        this->customPlot_Ascan->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
        //
        this->customPlot_Ascan->replot();
    }
}

