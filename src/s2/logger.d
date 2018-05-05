module s2.logger;

import s2.util.log.logger;

shared auto logger = new Logger!(LogLevel.ERROR)();
