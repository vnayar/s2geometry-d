module s2.util.log.logger;

import std.stdio;
import core.vararg;

/**
 * The log level is used to control which log statements should result in output.
 * Log levels are presumed to be inclusive, that is, ERROR implies FATAL and WARN
 * implies ERROR, FATAL, etc.
 */
enum LogLevel {
  /// No events are recorded.
  OFF,
  // Very severe error events that will presumably cause the application to abort.
  FATAL,
  /// Error events that may still allow the application to continue running.
  ERROR,
  /// Events that indicate a potentially harmful situation.
  WARN,
  /// Informational messages that highlight the progress of the application at a coarse level.
  INFO,
  /// Fine grain informational events useful for debugging.
  DEBUG,
  /// Finer grain information than DEBUG.
  TRACE
}

class Logger(LogLevel LogLevelV = LogLevel.DEBUG) {
public:

  // TODO: Add the optional ability to write to a file.
  this() {
  }

  static void log(LogLevel logLevel = LogLevel.DEBUG, T...)(lazy T args) {
    static if (logLevel <= LogLevelV) {
      writeln(args);
    }
  }

  static void logFatal(T...)(lazy T args) {
    log!(LogLevel.FATAL)(args);
  }

  static void logError(T...)(lazy T args) {
    log!(LogLevel.ERROR)(args);
  }

  static void logWarn(T...)(lazy T args) {
    log!(LogLevel.WARN)(args);
  }

  static void logDebug(T...)(lazy T args) {
    log!(LogLevel.DEBUG)(args);
  }

  static void logInfo(T...)(lazy T args) {
    log!(LogLevel.INFO)(args);
  }

  static void logTrace(T...)(lazy T args) {
    log!(LogLevel.TRACE)(args);
  }

  static void logf(LogLevel logLevel = LogLevel.DEBUG, T...)(string fmt, lazy T args) {
    static if (logLevel <= LogLevelV) {
      writefln(fmt, args);
    }
  }

  static void logfFatal(T...)(lazy T args) {
    logf!(LogLevel.FATAL)(args);
  }

  static void logfError(T...)(lazy T args) {
    logf!(LogLevel.ERROR)(args);
  }

  static void logfWarn(T...)(lazy T args) {
    logf!(LogLevel.WARN)(args);
  }

  static void logfInfo(T...)(lazy T args) {
    logf!(LogLevel.INFO)(args);
  }

  static void logfDebug(T...)(lazy T args) {
    logf!(LogLevel.DEBUG)(args);
  }

  static void logfTrace(T...)(lazy T args) {
    logf!(LogLevel.TRACE)(args);
  }

}

unittest {
  auto logger = new Logger!(LogLevel.DEBUG)();

  logger.log(typeid(logger));
  int i = 0;
  logger.logError("Hello error.", i++);
  assert(i == 1);
  logger.logDebug("Hello debug.", i++);
  assert(i == 2);
  logger.logTrace("Hello trace.", i++);
  assert(i == 2);
}

unittest {
  auto logger = new Logger!(LogLevel.ERROR)();

  logger.log(typeid(logger));
  int i = 0;
  logger.logFatal("Hello fatal.", i++);
  assert(i == 1);
  logger.logError("Hello error.", i++);
  assert(i == 2);
  logger.logDebug("Hello debug.", i++);
  assert(i == 2);
}

unittest {
  auto logger = new Logger!(LogLevel.OFF)();

  logger.log(typeid(logger));
  int i = 0;
  logger.logFatal("Hello fatal.", i++);
  assert(i == 0);
  logger.logError("Hello error.", i++);
  assert(i == 0);
  logger.logDebug("Hello debug.", i++);
  assert(i == 0);
}
