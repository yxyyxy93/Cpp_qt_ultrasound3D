#include <QVector>
#include <QWidget>
#include <QVBoxLayout>
#include <QCustomPlot>

class SliceWidget : public QWidget {
    Q_OBJECT
public:
    explicit SliceWidget(QVector<QVector<QVector<double>>> data, QWidget *parent = nullptr) : QWidget(parent), data_(data) {
        // Create a QCustomPlot widget and add it to the layout
        plot_ = new QCustomPlot();
        QVBoxLayout *layout = new QVBoxLayout(this);
        layout->addWidget(plot_);

        // Create an image plot and add it to the plot widget
        image_ = new QCPItemPixmap(plot_);
        plot_->addItem(image_);

        // Set the position and size of the image plot
        image_->setScaled(true);
        image_->setPixmap(QPixmap(data_.size(), data_[0].size()));
        image_->topLeft->setType(QCPItemPosition::ptAxisRectRatio);
        image_->topLeft->setCoords(0, 1);
        image_->bottomRight->setType(QCPItemPosition::ptAxisRectRatio);
        image_->bottomRight->setCoords(1, 0);

        // Connect the signal for updating the slice to the slot for updating the plot
        connect(this, &SliceWidget::sliceChanged, this, &SliceWidget::updatePlot);

        // Initialize the plot with the middle slice in the z-direction
        sliceIndex_ = data_[0][0].size()/2;
        emit sliceChanged();
    }

signals:
    void sliceChanged();

public slots:
    void setSliceIndex(int sliceIndex) {
        sliceIndex_ = sliceIndex;
        emit sliceChanged();
    }

private:
    void updatePlot() {
        // Get the slice you want to visualize
        int n1 = data_.size();
        int n2 = data_[0].size();
        QVector<QVector<double>> slice(n1, QVector<double>(n2));
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                slice[i][j] = data_[i][j][sliceIndex_];
            }
        }

        // Set the image data to the slice
        QImage image(n1, n2, QImage::Format_RGB32);
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                double value = slice[i][j];
                image.setPixel(i, j, qRgb(value, value, value));
            }
        }
        image = image.mirrored(false, true); // Flip the image vertically
        image = image.scaled(n1, n2); // Scale the image to match the size of the pixmap
        image = image.convertToFormat(QImage::Format_ARGB32); // Convert the image format to match the pixmap
        image.setDevicePixelRatio(plot_->devicePixelRatio()); // Set the device pixel ratio to match the plot widget
        image = image.copy(); // Make a copy of the image to prevent issues with shared data

        image_->setPixmap(QPixmap::fromImage(image)); // Set the pixmap to the image plot
        plot_->replot(); // Update the plot
    }
}
