#ifndef AMREX_ANSI_ESC_CODE_H_
#define AMREX_ANSI_ESC_CODE_H_

// https://en.wikipedia.org/wiki/ANSI_escape_code

namespace amrex {

constexpr char ResetDisplay[] = "\033[0m";

namespace Font {
    constexpr char Bold        [] = "\033[1m";
    constexpr char Faint       [] = "\033[2m";
    constexpr char Italic      [] = "\033[3m";
    constexpr char Underline   [] = "\033[4m";
    constexpr char SlowBlink   [] = "\033[5m";
    constexpr char RapidBlink  [] = "\033[6m";
}

namespace FGColor { // Foreground colors
    constexpr char Black         [] = "\033[30m";
    constexpr char Red           [] = "\033[31m";
    constexpr char Green         [] = "\033[32m";
    constexpr char Yellow        [] = "\033[33m";
    constexpr char Blue          [] = "\033[34m";
    constexpr char Magenta       [] = "\033[35m";
    constexpr char Cyan          [] = "\033[36m";
    constexpr char White         [] = "\033[37m";

    constexpr char Gray          [] = "\033[90m";
    constexpr char BrightBlack   [] = "\033[90m";
    constexpr char BrightRed     [] = "\033[91m";
    constexpr char BrightGreen   [] = "\033[92m";
    constexpr char BrightYellow  [] = "\033[93m";
    constexpr char BrightBLUE    [] = "\033[94m";
    constexpr char BrightMagenta [] = "\033[95m";
    constexpr char BrightCyan    [] = "\033[96m";
    constexpr char BrightWhite   [] = "\033[97m";
}

namespace BGColor { // Background colors
    constexpr char Black          [] = "\033[40m";
    constexpr char Red            [] = "\033[41m";
    constexpr char Green          [] = "\033[42m";
    constexpr char Yellow         [] = "\033[43m";
    constexpr char Blue           [] = "\033[44m";
    constexpr char Magenta        [] = "\033[45m";
    constexpr char Cyan           [] = "\033[46m";
    constexpr char White          [] = "\033[47m";

    constexpr char Gray           [] = "\033[100m";
    constexpr char BrightBlack    [] = "\033[100m";
    constexpr char BrightRed      [] = "\033[101m";
    constexpr char BrightGreen    [] = "\033[102m";
    constexpr char BrightYellow   [] = "\033[103m";
    constexpr char BrightBlue     [] = "\033[104m";
    constexpr char BrightMagenta  [] = "\033[105m";
    constexpr char BrightCyan     [] = "\033[106m";
    constexpr char BrightWhite    [] = "\033[107m";
}

}

#endif
