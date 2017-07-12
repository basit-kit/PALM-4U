//----------------------------------------------------------------------------//
// Current revisions:
// -----------------
// 
//
// Former revisions:
// -----------------
// $Id: mainwindow.h 1612 2015-07-07 12:25:21Z maronga $
//
// 1611 2015-07-07 12:23:22Z maronga
// Added slot start_watchdog 
//
// 793 2011-12-12 15:15:24Z maronga
// Initial revision
//
// Description:
// ------------
// UI mainwindow header file
//----------------------------------------------------------------------------//

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>


namespace Ui {
    class MainWindow;
}



class MainWindow : public QMainWindow
{
    Q_OBJECT


public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    int change_commandline(QString id,QString fwstring);
    int delete_commandline(QString id);
    int activate_flag(QString id);
    int deactivate_flag(QString id);
    int setup_gui(QString mrunline);
    int recent_jobs(int noj);


private slots:
    int startmrun();
    int enable_advanced();
    int enable_coupled();
    int choosejob();
    int choosejob_list();
    int change_rc_list();
    int check_flags();
    int change_lineinput();
    int reset_window();
    int save_to_file();
    int save_default();
    int start_watchdog();
    int open_from_file();
    int open_last();
    int help();
    int about_gui();
};

#endif // MAINWINDOW_H


